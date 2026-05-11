function fcdm = icdm_compute_functional_cdm(subj, K, opts)
% ICDM_COMPUTE_FUNCTIONAL_CDM  Compute voxelwise functional CDM from rs-fMRI.
%
%   fcdm = icdm_compute_functional_cdm(subj, K, opts)
%
%   For each white-matter voxel, computes a functional connectivity
%   distribution to K cortical ROIs from resting-state fMRI, yielding a
%   compositional profile analogous to the structural iCDM.
%
%   Pipeline:
%     1. Load rs-fMRI 4-D volume (native space)
%     2. Load cortical parcellation (same DK/VEP atlas as structural CDM)
%     3. Extract WM voxel BOLD time series at idx_native locations
%     4. Extract mean cortical ROI time series (K regions)
%     5. Compute Pearson correlation between each WM voxel and each ROI
%     6. Rectify (max(r,0)) and normalise to compositional profile
%
%   INPUT
%     subj : subject struct (from icdm_compose_subject)
%       Required fields:
%         .datafile    : cached MAT with idx_native, dim_native
%         .fmri_4d     : path to rs-fMRI 4-D NIfTI (native space,
%                        coregistered to DWI/T1)
%         .parcellation: path to cortical parcellation NIfTI (native space,
%                        integer labels 1..K)
%       Optional fields:
%         .Lesion      : [X x Y x Z] lesion mask (excluded from ROI means)
%     K    : number of cortical ROIs (e.g. 68 for DK, 162 for VEP)
%     opts : options struct
%       .fcdm.confounds    : path to confound regressors file (optional)
%       .fcdm.bandpass     : [fmin fmax] in Hz (default [0.01 0.1])
%       .fcdm.tr           : repetition time in seconds (required)
%       .fcdm.min_corr     : minimum correlation to retain (default 0)
%       .fcdm.smooth_fwhm  : spatial smoothing FWHM in mm (default 0, none)
%       .idx_regions       : region index selection (default 1:K)
%
%   OUTPUT
%     fcdm : struct with fields
%       .pi_func     : [Nnat x K] functional CDM (compositional, rows sum to 1)
%       .corr_raw    : [Nnat x K] raw Pearson correlations (before rectification)
%       .idx_native  : [Nnat x 1] native voxel indices
%       .dim_native  : [1 x 3] native volume dimensions
%       .roi_tc      : [T x K] mean ROI time courses (for inspection)
%
%   See also ICDM_SUBJECT_VB, ICDM_VALIDATE_STRUCTURE_FUNCTION

fprintf('[fCDM] Subject %s\n', subj.id);

% ---- defaults ----
fcdm_opts = getfield_default(opts, 'fcdm', struct());
TR         = fcdm_opts.tr;           % must be provided
bandpass   = getfield_default(fcdm_opts, 'bandpass', [0.01 0.1]);
min_corr   = getfield_default(fcdm_opts, 'min_corr', 0);
smooth_fwhm= getfield_default(fcdm_opts, 'smooth_fwhm', 0);
idx_regions= getfield_default(opts, 'idx_regions', 1:K);
K_sel      = numel(idx_regions);

% ---- 1. Load cached native-space info ----
load(subj.datafile, 'idx_native', 'dim_native');
Nnat = numel(idx_native);
[Xn, Yn, Zn] = deal(dim_native(1), dim_native(2), dim_native(3));

% ---- 2. Load rs-fMRI ----
fprintf('  Loading rs-fMRI: %s\n', subj.fmri_4d);
Vf = spm_vol(subj.fmri_4d);
T  = numel(Vf);
fprintf('  %d volumes, TR=%.2f s\n', T, TR);

% Read all volumes into [X Y Z T]
fmri_4d = spm_read_vols(Vf);

% Optional spatial smoothing
if smooth_fwhm > 0
    fprintf('  Smoothing fMRI (FWHM=%.1f mm)...\n', smooth_fwhm);
    for t = 1:T
        spm_smooth(fmri_4d(:,:,:,t), fmri_4d(:,:,:,t), smooth_fwhm);
    end
end

% ---- 3. Load parcellation ----
fprintf('  Loading parcellation: %s\n', subj.parcellation);
Vp = spm_vol(subj.parcellation);
parc_vol = round(spm_read_vols(Vp));  % integer labels

% Lesion mask (optional): exclude lesion voxels from ROI means
if isfield(subj, 'Lesion') && ~isempty(subj.Lesion)
    lesion_mask = subj.Lesion > 0.5;
else
    lesion_mask = false(Xn, Yn, Zn);
end

% ---- 4. Extract ROI mean time series ----
fprintf('  Extracting %d ROI time series...\n', K_sel);
roi_tc = zeros(T, K_sel, 'single');
roi_nvox = zeros(1, K_sel);

for k = 1:K_sel
    roi_idx = find(parc_vol(:) == idx_regions(k) & ~lesion_mask(:));
    roi_nvox(k) = numel(roi_idx);
    if roi_nvox(k) < 5
        warning('ROI %d has only %d voxels (< 5)', idx_regions(k), roi_nvox(k));
    end
    if roi_nvox(k) > 0
        roi_ts = reshape(fmri_4d, [], T);
        roi_tc(:,k) = single(mean(roi_ts(roi_idx, :), 1)');
    end
end

% ---- 5. Extract WM voxel time series ----
fprintf('  Extracting %d WM voxel time series...\n', Nnat);
wm_tc = zeros(T, Nnat, 'single');
fmri_2d = reshape(fmri_4d, [], T);
wm_tc = single(fmri_2d(idx_native, :)');   % [T x Nnat]
clear fmri_4d fmri_2d;

% ---- 6. Confound regression (optional) ----
if isfield(fcdm_opts, 'confounds') && ~isempty(fcdm_opts.confounds)
    fprintf('  Regressing confounds...\n');
    R = load(fcdm_opts.confounds);
    if isstruct(R), R = struct2array(R); end
    R = [ones(T,1), R];   % add intercept
    % Project out confounds from both WM and ROI time series
    proj = eye(T) - R * (R \ eye(T));
    wm_tc  = proj * wm_tc;
    roi_tc = proj * roi_tc;
end

% ---- 7. Bandpass filter ----
if ~isempty(bandpass) && all(bandpass > 0)
    fprintf('  Bandpass filtering [%.3f - %.3f] Hz...\n', bandpass(1), bandpass(2));
    wm_tc  = bandpass_filter(wm_tc, TR, bandpass(1), bandpass(2));
    roi_tc = bandpass_filter(roi_tc, TR, bandpass(1), bandpass(2));
end

% ---- 8. Compute Pearson correlations ----
fprintf('  Computing voxelwise correlations...\n');
% Z-score time series
wm_z  = zscore(wm_tc, 0, 1);   % [T x Nnat]
roi_z = zscore(roi_tc, 0, 1);   % [T x K_sel]

% Correlation matrix: [Nnat x K_sel]
corr_raw = single((wm_z' * roi_z) / (T - 1));
corr_raw(~isfinite(corr_raw)) = 0;

% ---- 9. Rectify and normalise to compositional profile ----
corr_rect = max(corr_raw, min_corr);
row_sum = sum(corr_rect, 2);
row_sum(row_sum < eps) = 1;   % avoid division by zero
pi_func = corr_rect ./ row_sum;

% Voxels with all-zero correlations: set to uniform
uniform = single(ones(1, K_sel) / K_sel);
zero_rows = row_sum < eps;
pi_func(zero_rows, :) = repmat(uniform, sum(zero_rows), 1);

fprintf('  fCDM computed: %d voxels x %d ROIs\n', Nnat, K_sel);
fprintf('  Mean max-correlation: %.3f\n', mean(max(corr_raw, [], 2)));
fprintf('  Voxels with all-zero functional profile: %d (%.1f%%)\n', ...
    sum(zero_rows), 100*sum(zero_rows)/Nnat);

% ---- output ----
fcdm.pi_func    = pi_func;
fcdm.corr_raw   = corr_raw;
fcdm.idx_native = idx_native;
fcdm.dim_native = dim_native;
fcdm.roi_tc     = roi_tc;
fcdm.roi_nvox   = roi_nvox;

end


%% ---- LOCAL: simple bandpass via FFT ----
function Y = bandpass_filter(Y, TR, fmin, fmax)
% Bandpass filter columns of Y using FFT.
%   Y : [T x N]  time series (columns)
%   TR : repetition time (seconds)
%   fmin, fmax : passband edges (Hz)

[T, ~] = size(Y);
freq = (0:T-1)' / (T * TR);   % frequency axis

% Two-sided mask
mask = (freq >= fmin & freq <= fmax) | (freq >= (1/TR - fmax) & freq <= (1/TR - fmin));
mask(1) = 0;   % remove DC

Y_fft = fft(Y, [], 1);
Y_fft(~mask, :) = 0;
Y = real(ifft(Y_fft, [], 1));
end
