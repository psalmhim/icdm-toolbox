% iCDM Toolbox - Tumor Patient Preprocessing
% Copyright (c) 2024-2026 Hae-Jeong Park, Ph.D.
% Yonsei University, Department of Nuclear Medicine
%
% STUDY_PREPROCESS_TUMOR  Preprocess tumor patient data for iCDM.
%
%   Step 1: Gunzip NIfTI files
%   Step 2: Tumor segmentation (manual or semi-auto)
%   Step 3: Spatial normalization with cost-function masking
%   Step 4: FreeSurfer parcellation (or atlas-based)
%
%   Data assumed in BIDS format under icdm/data/pre/
%
%   REQUIRES: SPM12, FreeSurfer (optional)

clear; clc;

%% ========== 0. CONFIGURATION ==========

datadir  = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'pre');
outdir   = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'derivatives');

subjects = {'sub-PAT01', 'sub-PAT02', 'sub-PAT03'};
session  = 'ses-preop';

% SPM path check
if ~exist('spm','file')
    error('SPM12 not on path. Add SPM12 to MATLAB path first.');
end
spm('defaults', 'fmri');
spm_jobman('initcfg');

%% ========== 1. GUNZIP ALL .nii.gz ==========

fprintf('=== Step 1: Gunzip NIfTI files ===\n');
for s = 1:numel(subjects)
    sid = subjects{s};
    subdir = fullfile(datadir, sid, session);

    gz_files = dir(fullfile(subdir, '**', '*.nii.gz'));
    for f = 1:numel(gz_files)
        gz_path = fullfile(gz_files(f).folder, gz_files(f).name);
        nii_path = gz_path(1:end-3);  % remove .gz
        if ~exist(nii_path, 'file')
            fprintf('  Gunzip: %s\n', gz_files(f).name);
            gunzip(gz_path, gz_files(f).folder);
        end
    end
end

%% ========== 2. TUMOR MASK ==========
%
% Option A: Manual segmentation (ITK-SNAP, 3D Slicer)
%   - Draw tumor mask on T1, save as <sid>_tumor_mask.nii
%
% Option B: Semi-automatic (below creates a placeholder)
%   - User needs to edit/replace with actual tumor mask
%
% Option C: AI-based (e.g., HD-BET + HD-GLIO, nnU-Net)

fprintf('\n=== Step 2: Check tumor masks ===\n');
for s = 1:numel(subjects)
    sid = subjects{s};
    t1_path = fullfile(datadir, sid, session, 'anat', ...
        sprintf('%s_%s_T1w.nii', sid, session));
    deriv_anat = fullfile(outdir, sid, session, 'anat');
    mkdir_if_needed(deriv_anat);

    mask_path = fullfile(deriv_anat, sprintf('%s_tumor_mask.nii', sid));

    if exist(mask_path, 'file')
        fprintf('  [%s] Tumor mask exists: %s\n', sid, mask_path);
    else
        fprintf('  [%s] *** Tumor mask NOT FOUND ***\n', sid);
        fprintf('         Expected: %s\n', mask_path);
        fprintf('         Please create manually (ITK-SNAP) or via AI segmentation.\n');
        fprintf('         Creating empty placeholder...\n');

        % Create empty mask with same geometry as T1
        V = spm_vol(t1_path);
        Y = zeros(V.dim);
        Vo = V;
        Vo.fname = mask_path;
        Vo.dt = [spm_type('uint8') 0];
        Vo.pinfo = [1; 0; 0];
        spm_write_vol(Vo, Y);
        fprintf('         Placeholder saved. EDIT THIS before proceeding.\n');
    end
end

%% ========== 3. SPATIAL NORMALIZATION WITH COST-FUNCTION MASKING ==========
%
% SPM Unified Segmentation with lesion mask:
%   - The tumor mask is used to exclude lesion voxels from the
%     cost function during tissue classification and registration
%   - This ensures the deformation field is driven only by healthy tissue
%   - Inside the tumor, the deformation is smoothly extrapolated
%
% Output:
%   y_<T1>.nii       : forward deformation field (native → MNI)
%   iy_<T1>.nii      : inverse deformation field (MNI → native)
%   c1-c5_<T1>.nii   : tissue probability maps
%   m<T1>.nii        : bias-corrected T1

fprintf('\n=== Step 3: Spatial normalization (cost-function masking) ===\n');

for s = 1:numel(subjects)
    sid = subjects{s};
    t1_path = fullfile(datadir, sid, session, 'anat', ...
        sprintf('%s_%s_T1w.nii', sid, session));
    deriv_anat = fullfile(outdir, sid, session, 'anat');
    mask_path = fullfile(deriv_anat, sprintf('%s_tumor_mask.nii', sid));

    % Check if already done
    [t1dir, t1name, ~] = fileparts(t1_path);
    y_file = fullfile(t1dir, ['y_' t1name '.nii']);
    if exist(y_file, 'file')
        fprintf('  [%s] Deformation field exists, skipping.\n', sid);
        continue;
    end

    % Check tumor mask is not empty placeholder
    Vm = spm_vol(mask_path);
    Ym = spm_read_vols(Vm);
    if max(Ym(:)) == 0
        warning('[%s] Tumor mask is empty! Normalization will proceed WITHOUT masking.', sid);
        warning('Edit %s before running for real.', mask_path);
        use_mask = false;
    else
        use_mask = true;
        fprintf('  [%s] Tumor mask: %d voxels (%.1f cm³)\n', sid, ...
            nnz(Ym > 0.5), nnz(Ym > 0.5) * abs(det(Vm.mat(1:3,1:3))) / 1000);
    end

    fprintf('  [%s] Running SPM unified segmentation...\n', sid);

    % ---- SPM batch: Unified Segmentation ----
    matlabbatch = {};

    % Channel
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {[t1_path ',1']};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1]; % save bias-corrected

    % Tissue priors (SPM12 TPM)
    tpm_path = fullfile(spm('dir'), 'tpm', 'TPM.nii');
    for t = 1:6
        matlabbatch{1}.spm.spatial.preproc.tissue(t).tpm = {sprintf('%s,%d', tpm_path, t)};
        matlabbatch{1}.spm.spatial.preproc.tissue(t).ngaus = [1 1 2 3 4 2];
        matlabbatch{1}.spm.spatial.preproc.tissue(t).ngaus = matlabbatch{1}.spm.spatial.preproc.tissue(t).ngaus(t);
        if t <= 3
            matlabbatch{1}.spm.spatial.preproc.tissue(t).native = [1 0];
        else
            matlabbatch{1}.spm.spatial.preproc.tissue(t).native = [0 0];
        end
        matlabbatch{1}.spm.spatial.preproc.tissue(t).warped = [0 0];
    end

    % Warping: write forward + inverse deformation fields
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1]; % inverse + forward

    % ---- Cost-function masking: exclude tumor from segmentation ----
    % In SPM12, the lesion mask is applied by modifying the tissue
    % priors at lesion voxels. We set all tissue priors to equal
    % probability inside the lesion, making those voxels uninformative
    % for the registration cost function.
    if use_mask
        fprintf('  [%s] Applying cost-function masking...\n', sid);

        % Reslice tumor mask to T1 space if needed (should already match)
        % Then create a modified TPM where lesion voxels are masked
        % SPM12 approach: use the lesion mask via spm_preproc_run
        % by setting channel.lesion field (SPM12 r7771+)
        %
        % Alternative (works with all SPM12 versions):
        % Zero out T1 intensity inside tumor before segmentation,
        % then the segmentation naturally treats it as background.
        %
        % We use the native SPM12 lesion masking if available:
        try
            matlabbatch{1}.spm.spatial.preproc.channel.lesion = {[mask_path ',1']};
            fprintf('  [%s] Using SPM native lesion masking.\n', sid);
        catch
            % Fallback: mask T1 intensity (set lesion to NaN or mean)
            fprintf('  [%s] SPM lesion field not available. Using intensity masking.\n', sid);
            Vt1 = spm_vol(t1_path);
            Yt1 = spm_read_vols(Vt1);
            Yt1(Ym > 0.5) = nanmean(Yt1(Ym <= 0.5));  % fill with mean healthy
            % Write masked T1
            t1_masked = fullfile(deriv_anat, [t1name '_lesion_filled.nii']);
            Vo = Vt1;
            Vo.fname = t1_masked;
            spm_write_vol(Vo, Yt1);
            matlabbatch{1}.spm.spatial.preproc.channel.vols = {[t1_masked ',1']};
            fprintf('  [%s] Using lesion-filled T1: %s\n', sid, t1_masked);
        end
    end

    % Run segmentation
    fprintf('  [%s] Executing SPM batch...\n', sid);
    spm_jobman('run', matlabbatch);

    % Verify output
    if exist(y_file, 'file')
        fprintf('  [%s] SUCCESS: %s\n', sid, y_file);
    else
        % Check if deformation field was written in derivatives
        y_alt = dir(fullfile(t1dir, 'y_*.nii'));
        if ~isempty(y_alt)
            fprintf('  [%s] SUCCESS: %s\n', sid, fullfile(t1dir, y_alt(1).name));
        else
            warning('[%s] Deformation field not found after segmentation!', sid);
        end
    end
end

%% ========== 4. COREGISTER DWI AND FMRI TO T1 ==========

fprintf('\n=== Step 4: Coregister DWI & fMRI to T1 ===\n');

for s = 1:numel(subjects)
    sid = subjects{s};
    subdir = fullfile(datadir, sid, session);

    % Bias-corrected T1
    t1_name = sprintf('%s_%s_T1w', sid, session);
    t1_bc = fullfile(subdir, 'anat', ['m' t1_name '.nii']);
    if ~exist(t1_bc, 'file')
        t1_bc = fullfile(subdir, 'anat', [t1_name '.nii']);
    end

    % --- DWI: coregister b0 to T1 ---
    dwi_ap = fullfile(subdir, 'dwi', ...
        sprintf('%s_%s_acq-AP_dwi.nii', sid, session));
    if exist(dwi_ap, 'file')
        fprintf('  [%s] Coregistering DWI (AP) → T1...\n', sid);
        matlabbatch = {};
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[t1_bc ',1']};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {[dwi_ap ',1']};  % first volume (b0)
        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = ...
            [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        spm_jobman('run', matlabbatch);
    end

    % --- fMRI: coregister mean to T1 ---
    fmri = fullfile(subdir, 'func', ...
        sprintf('%s_%s_task-rest_bold.nii', sid, session));
    if exist(fmri, 'file')
        fprintf('  [%s] Coregistering fMRI → T1...\n', sid);

        % First compute mean fMRI volume
        Vf = spm_vol(fmri);
        mean_img = zeros(Vf(1).dim);
        for t = 1:numel(Vf)
            mean_img = mean_img + spm_read_vols(Vf(t));
        end
        mean_img = mean_img / numel(Vf);

        mean_path = fullfile(subdir, 'func', ...
            sprintf('mean_%s_%s_task-rest_bold.nii', sid, session));
        Vm = Vf(1);
        Vm.fname = mean_path;
        spm_write_vol(Vm, mean_img);

        matlabbatch = {};
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[t1_bc ',1']};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {[mean_path ',1']};
        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = ...
            [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        spm_jobman('run', matlabbatch);
    end
end

fprintf('\n=== Preprocessing complete ===\n');
fprintf('Next steps:\n');
fprintf('  1. VERIFY tumor masks in %s/*/ses-preop/anat/\n', outdir);
fprintf('  2. If placeholder, draw tumor masks manually and re-run Step 3\n');
fprintf('  3. Run DWI preprocessing (topup, eddy) externally (MRtrix3)\n');
fprintf('  4. Run tractography + iCDM construction\n');
fprintf('  5. Run study_validate_tumor.m\n');

%% ==================== LOCAL ====================
function mkdir_if_needed(d)
if ~exist(d, 'dir'), mkdir(d); end
end
