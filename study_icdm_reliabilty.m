% iCDM Toolbox - Individualized Compositional Diffusion Microstructure
% Copyright (c) 2024-2026 Hae-Jeong Park, Ph.D.
% Yonsei University, Department of Nuclear Medicine
%
% This software is part of the iCDM framework described in:
%   Park et al., "Individualized Connection Distribution Mapping:
%   A Hierarchical Bayesian Framework for Voxelwise Compositional
%   Connectivity Inference from Diffusion MRI" (submitted).
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

% STUDY_ICDM_RELIABILTY  Script for split-half reliability analysis of ICDM data.
%
%   Processes multiple HBN subjects through dodti_icdm_process_reliability
%   to generate split-half fold files, computes voxelwise KL and JS
%   divergence between folds (icdm_cdm_kl_js), warps results to MNI
%   template space, and loads the group mask for further analysis.
%
%   See also ICDM_CDM_KL_JS, GENERIC_WARP_4D
% Author: Hae-Jeong Park, Ph.D.

%subjects(12)
if 0
subjpath='/remotenas2/projects/hbn_connectome/vbc_tck3_10M/OUT/SC/HBN/NDARDK983BDA/';
icdmfullfile='/remotenas2/projects/hbn_connectome/vbc_tck3_10M/OUT/SC/HBN/NDARDK983BDA/WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm.nii';
vfull=spm_vol(icdmfullfile);
tckfile='WBT_10M_ctx.tck'; skipmode=1; parcidx=[1:84];
lab='DesikanKilliany_68Parcels_to_fs_t1_native+subctx.nii';
labfiles=dodti_icdm_process_reliability(subjpath,tckfile,lab, skipmode,0,parcidx);

%subjects(10)
subjpath='/remotenas2/projects/hbn_connectome/vbc_tck3_10M/OUT/SC/HBN/NDARDH086ZKK/';
icdmfullfile='/remotenas2/projects/hbn_connectome/vbc_tck3_10M/OUT/SC/HBN/NDARDH086ZKK/WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm.nii';
vfull=spm_vol(icdmfullfile);
tckfile='WBT_10M_ctx.tck'; skipmode=1; parcidx=[1:84];
lab='DesikanKilliany_68Parcels_to_fs_t1_native+subctx.nii';
labfiles=dodti_icdm_process_reliability(subjpath,tckfile,lab, skipmode,0,parcidx);


subjpath='/remotenas2/projects/hbn_connectome/vbc_tck3_10M/OUT/SC/HBN/NDARAM277WZT/';
icdmfullfile='/remotenas2/projects/hbn_connectome/vbc_tck3_10M/OUT/SC/HBN/NDARAM277WZT/WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm.nii';
vfull=spm_vol(icdmfullfile);
tckfile='WBT_10M_ctx.tck'; skipmode=1; parcidx=[1:84];
lab='DesikanKilliany_68Parcels_to_fs_t1_native+subctx.nii';
labfiles=dodti_icdm_process_reliability(subjpath,tckfile,lab, skipmode,0,parcidx);

subjpath='/remotenas2/projects/hbn_connectome/vbc_tck3_10M/OUT/SC/HBN/NDARBF998MBA/';
icdmfullfile='/remotenas2/projects/hbn_connectome/vbc_tck3_10M/OUT/SC/HBN/NDARBF998MBA/WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm.nii';
vfull=spm_vol(icdmfullfile);
tckfile='WBT_10M_ctx.tck'; skipmode=1; parcidx=[1:84];
lab='DesikanKilliany_68Parcels_to_fs_t1_native+subctx.nii';
labfiles=dodti_icdm_process_reliability(subjpath,tckfile,lab, skipmode,0,parcidx);

end

%% ===== USER CONFIGURATION (edit before running) =====
% mpath = '/path/to/your/ICDM_desikan68';
error('Edit the USER CONFIGURATION section above before running this script.');

mpath='/Volumes/eHDD/DATA/ICDM_desikan68';

fold1_file = 'WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm_fold1.nii';
fold2_file = 'WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm_fold2.nii';
mask_r_file  = 'rbrain_277.nii';


[outKL, outJS] = icdm_cdm_kl_js(fullfile(mpath,fold1_file),fullfile(mpath,fold2_file),[],fullfile(mpath,'cdm_277'));

[outKL1, outJS1] = icdm_cdm_kl_js(fullfile(mpath,'P_age8.nii'),fullfile(mpath,'P_age12.nii'),[],fullfile(mpath,'cdm_age'));
reg_image(outKL1,outJS1)
reg_image(outKL,outJS)

mpath='../simpy/data/';
def_field = fullfile(mpath,'dartel','y_fs_t1.nii');
templatefile=fullfile(mpath,'dartel','HBN_Temp_6.nii'); Vtemplate=spm_vol(templatefile);
icdmfile1=fullfile(mpath,'WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm_fold1.nii');
Vicdm1=spm_vol(icdmfile1);
icdmfile2=fullfile(mpath,'WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm_fold2.nii');
Vicdm2=spm_vol(icdmfile2);


filenames={'cdm_277_JS.nii', 'cdm_277_KL12.nii','cdm_277_Nv.nii','cdm_age_JS.nii', 'cdm_age_KL12.nii','WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm_fold1.nii', 'WBT_10M_ctx_DesikanKilliany_68Parcels_to_fs_t1_native+subctx_to_dwi15_icdm_fold2.nii'};
for i=1:length(filenames)
v=spm_vol(fullfile(mpath,filenames{i}));
Yin=spm_read_vols(v);
Yout = generic_warp_4d(Yin, Vicdm1(1), def_field, templatefile, +1, 1);
vo=Vtemplate(1);
vo.fname=fullfile(mpath,['w' filenames{i}]);
vo.dt=v.dt;
vol_write_vols(vo,Yout);
end


load(fullfile(mpath,'group_mask.mat'))