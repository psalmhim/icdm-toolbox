% iCDM Toolbox - Individualized Compositional Diffusion Microstructure
% Copyright (c) 2024-2026 Hae-Jeong Park, Ph.D.
% Yonsei University, Department of Nuclear Medicine
%
% This software is part of the iCDM framework described in:
%   Park et al., "Individualized Connection Distribution Mapping:
%   A Hierarchical Bayesian Framework for Voxelwise Compositional
%   Connectivity Inference from Diffusion MRI", NeuroImage (2026).
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%   2. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
%   3. Neither the name of the copyright holder nor the names of its
%      contributors may be used to endorse or promote products derived
%      from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES.
%
% REQUIRES: SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)

% STUDY_SIM_RELIABILITY  Split-half stability evaluation of CDM using SPM.
%
%   Loads two split-half ICDM fold NIfTI files (volumes 1:68),
%   normalizes connectivity within the 68 cortical ROIs, reslices a
%   T1-derived brain mask to DWI space, and computes the voxelwise
%   Jensen-Shannon divergence as a measure of test-retest reliability.
%   Produces NIfTI outputs (JS map, Nv map) and a three-panel figure
%   showing sagittal, coronal, and axial slices at the maximum-JS
%   location, plus a log(Nv)-vs-log(JS) scatter plot.
%
%   See also ICDM_CDM_KL_JS, ICDM_KL_JS
% Author: Hae-Jeong Park, Ph.D.

%% ============================================================
% Split-half stability evaluation using SPM
% - spm_vol / spm_read_vols
% - volumes 1:68 only
% - normalize within 68 ROIs
% - reslice T1 mask -> DWI space
% - compute JS divergence
% - show max-difference sagittal / coronal / axial slices
%% ============================================================

clear; clc;

spm('Defaults','fmri');
spm_jobman('initcfg');

EPS = 1e-12;

%% ------------------------------------------------------------
% INPUT FILES
%% ------------------------------------------------------------
fold1_file = '/Volumes/eHDD/DATA/epilepsy/NDARDK983BDA/WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm_fold1.nii';

fold2_file = '/Volumes/eHDD/DATA/epilepsy/NDARDK983BDA/WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm_fold2.nii';

mask_r_file  = '/Volumes/eHDD/DATA/epilepsy/NDARDK983BDA/rfs_t1_brain.nii';

out_js = 'split_half_JS_68.nii';
out_nv = 'split_half_Nv_68.nii';

%% ------------------------------------------------------------
% LOAD FOLD DATA (DWI SPACE)
%% ------------------------------------------------------------
V1 = spm_vol(fold1_file);
V2 = spm_vol(fold2_file);
ref = V1(1);
% restrict to volumes 1:68 (MATLAB indexing)
V1 = V1(1:68);
V2 = V2(1:68);

[C1, ~] = spm_read_vols(V1);
[C2, ~] = spm_read_vols(V2);

[X,Y,Z,K] = size(C1);
fprintf('[INFO] Using volumes 1–68 only | shape = (%d,%d,%d,%d)\n', X,Y,Z,K);

Vmask = spm_vol(mask_r_file);
mask  = spm_read_vols(Vmask) > 0;

fprintf('[INFO] Mask resliced to DWI space\n');

%% ------------------------------------------------------------
% STREAMLINE COUNT (Nv within 68 ROIs)
%% ------------------------------------------------------------
Nv = sum(C1,4) + sum(C2,4);
Nv(~mask) = NaN;

%% ------------------------------------------------------------
% JS DIVERGENCE COMPUTATION
%% ------------------------------------------------------------
JS = NaN(X,Y,Z);

for ix = 1:X
    for iy = 1:Y
        for iz = 1:Z

            if ~mask(ix,iy,iz)
                continue;
            end

            c1 = squeeze(C1(ix,iy,iz,:));
            c2 = squeeze(C2(ix,iy,iz,:));

            if sum(c1)==0 || sum(c2)==0
                continue;
            end

            % normalize within 68 ROIs
            p = c1 / max(sum(c1), EPS);
            q = c2 / max(sum(c2), EPS);

            p = max(p, EPS);
            q = max(q, EPS);
            m = 0.5*(p + q);

            kl_pm = sum(p .* log(p ./ m));
            kl_qm = sum(q .* log(q ./ m));

            JS(ix,iy,iz) = 0.5*(kl_pm + kl_qm);
        end
    end
end

fprintf('[INFO] JS divergence computed\n');

%% ------------------------------------------------------------
% SAVE NIFTI OUTPUTS
%% ------------------------------------------------------------
Vout = ref;


Vout.fname   = out_nv;
Vout.descrip = 'Split-half streamline count (68 ROIs)';
spm_write_vol(Vout, Nv);

Vout.fname   = out_js;
Vout.dt(1)=16;
Vout.descrip = 'Split-half JS divergence (68 ROIs)';
spm_write_vol(Vout, JS);


fprintf('[SAVED] %s\n', out_js);
fprintf('[SAVED] %s\n', out_nv);

%% ------------------------------------------------------------
% FIND SLICES WITH MAX DIFFERENCE
%% ------------------------------------------------------------
JS0 = JS;
JS0(isnan(JS0)) = 0;

[~, x_idx] = max(mean(JS0,[2 3]));  % sagittal
[~, y_idx] = max(mean(JS0,[1 3]));  % coronal
[~, z_idx] = max(mean(JS0,[1 2]));  % axial

fprintf('\n[MAX JS SLICES]\n');
fprintf('  Sagittal (x) : %d\n', x_idx);
fprintf('  Coronal  (y) : %d\n', y_idx);
fprintf('  Axial    (z) : %d\n', z_idx);

%% ------------------------------------------------------------
% VISUALIZATION
%% ------------------------------------------------------------
figure('Color','w','Position',[100 100 1200 400]);

subplot(1,3,1);
imagesc(squeeze(JS(x_idx,:,:))');
axis image off; colormap hot; colorbar;
title(sprintf('Sagittal (x=%d)', x_idx));

subplot(1,3,2);
imagesc(squeeze(JS(:,y_idx,:))');
axis image off; colormap hot; colorbar;
title(sprintf('Coronal (y=%d)', y_idx));

subplot(1,3,3);
imagesc(squeeze(JS(:,:,z_idx))');
axis image off; colormap hot; colorbar;
title(sprintf('Axial (z=%d)', z_idx));

sgtitle('Split-half JS divergence (SPM-based)');

aa=Nv+JS;
id=isnan(aa) & Nv>100;
nv=Nv(~id); js=JS(~id);

figure;
plot(log(nv+1),log(js),'.');
xlabel('log(Nv)');ylabel('log(JS)')