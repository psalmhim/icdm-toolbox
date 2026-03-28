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

% CONVERT  Utility script to convert per-subject MNI results to compact format.
%
%   Loads each subject MAT file in iter_001/, extracts y_ilr_mni and
%   kappa_mni at the group mask indices (idx_mni), and overwrites the file
%   with the reduced representation.  This reduces disk usage for large
%   cohorts where the full 4-D volumes are not needed.
%
%   See also ICDM_4D_TO_2D, ICDM_2D_TO_4D
% Author: Hae-Jeong Park, Ph.D.

%% ===== USER CONFIGURATION (edit before running) =====
% mpath = '/path/to/your/ICDM_desikan68/';
error('Edit the USER CONFIGURATION section above before running this script.');

mpath='~/data/ICDM_desikan68/';
load(fullfile(mpath,'group_mask.mat'));
files=mnet_file_lists([mpath '/iter_001/*.mat']);
for i=1:size(files,1)
    fn=deblank(files(i,:));
    a=load(fn);
    b.id=a.OUT.id;
    b.y_ilr_mni=icdm_4d_to_2d(a.OUT.y_ilr_mni,idx_mni);
    b.kappa_mni=a.OUT.kappa_mni(idx_mni);
    save(fn,'-struct','b','-v7.3');
end
