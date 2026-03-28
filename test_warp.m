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

% TEST_WARP  Test script for the iCDM trilinear warp operator system.
%
%   Loads test.mat containing subjects and opts, builds forward and
%   backward trilinear warp operators for one subject, and verifies
%   round-trip accuracy by comparing the prebuilt-operator warp against
%   the SPM-based spm_deformations warp.  Outputs are written as NIfTI
%   files for visual comparison via spm_check_registration.
%
%   See also ICDM_BUILD_WARP_TO_MNI, ICDM_BUILD_WARP_TO_NATIVE,
%            ICDM_WARP_TO_MNI, ICDM_WARP_TO_NATIVE, ICDM_WARP_4D
% Author: Hae-Jeong Park, Ph.D.

load test.mat

% Required fields
subj=subjects(100);
subj.dartel_flow = subj.dartel_flow;
subj.template    = subj.template;
subj.icdm_4d     = subj.icdm_4d;

% Group MNI grid info
grp.idx_mni = opts.idx_mni;  % linear index
grp.idx_mni = grp.idx_mni(:);
grp.dim_mni = opts.dim_mni;            % example
vt=spm_vol(subj.template);
Vn = spm_vol(subj.icdm_4d);
[p,f,e] = fileparts(subj.icdm_4d);
vidl=spm_vol(fullfile(p,[f '_idl' e]));
R=spm_read_vols(vidl);


fprintf('Building warp operators...\n');

subj.warp.to_native = icdm_build_warp_to_native( ...
    subj.icdm_4d, subj.dartel_flow,subj.template, opts);

subj.warp.to_mni = icdm_build_warp_to_mni( ...
    subj.icdm_4d, subj.dartel_flow, subj.template, opts);

Ywmni = icdm_warp_to_mni(R, subj.warp);
vo=vt(1);
vo.fname='wtestwarp.nii';  
spm_write_vol(vo, Ywmni);

Ywind = icdm_warp_to_native(Ywmni, subj.warp);
vo=Vn(1);
vo.fname='wtestwarp_ind.nii';
spm_write_vol(vo, Ywind);

%%
Yspm = icdm_warp_4d(R, ...
        Vn(1), subj.dartel_flow, subj.template, 1, 1);
vo=vt(1);
vo.fname='testwarp.nii';  
spm_write_vol(vo, Yspm);

Yind = icdm_warp_4d(Yspm, ...
        Vn(1), subj.dartel_flow, subj.template, -1, 1);

vo=Vn(1);
vo.fname='testwarp_ind.nii';
spm_write_vol(vo, Yind);

spm_check_registration('single_subj_T1.nii','testwarp.nii','wtestwarp.nii')
spm_check_registration('native.nii','testwarp_ind.nii','wtestwarp_ind.nii')