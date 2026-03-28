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

% STUDY_APPLY_TO_NEW  Apply a trained iCDM model to new subjects for individual prediction.
%
%   Loads the group mask and iteration results, then runs
%   icdm_predict_individual_cdm on a selected subject to produce
%   group-predicted (pi_grp), individual-posterior (pi_ind), and fused
%   (pi_fused) connectivity maps.
%
%   See also ICDM_PREDICT_INDIVIDUAL_CDM, ICDM_RESULT_TO_NII
% Author: Hae-Jeong Park, Ph.D.

cd(outdir);
load group_mask.mat
icdm_result_to_nii(outdir);
nid=90;
subj=subjects(nid);
sid  = subj.id;
subjfile = fullfile(outdir,'iter_003',[sid '_vb.mat']);
VBs = load(subjfile);
grp=load(fullfile(outdir,'group_icdm_iter_003.mat'));
VBc = icdm_predict_individual_cdm(subj, VBs, grp, opts);
