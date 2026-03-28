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

function subjects = icdm_build_subject_warps(subjects)
% ICDM_BUILD_SUBJECT_WARPS  Build or load trilinear warp operators for all subjects.
%
%   subjects = icdm_build_subject_warps(subjects)
%
%   Iterates over the subject struct array and ensures that each subject
%   has precomputed forward (native-to-MNI) and backward (MNI-to-native)
%   trilinear warp operators stored in its data file.  If the data file
%   already exists the warp is loaded; otherwise the ICDM 4-D volume is
%   read, a WM mask is computed, and both warp operators are built via
%   icdm_build_warp_to_native / icdm_build_warp_to_mni and saved.
%
%   Input / Output
%     subjects : struct array produced by icdm_compose_subject.
%                On return, each element has a valid datafile containing
%                the warp operators.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_COMPOSE_SUBJECT, ICDM_BUILD_WARP_TO_MNI,
%            ICDM_BUILD_WARP_TO_NATIVE

S = numel(subjects);
if S == 0
    error('No subjects given.');
end

for s=1:S
    if rem(s,10)==0
        fprintf('  check subject datafile %d/%d: %s\n', s, S, subjects(s).id);
    end
    subj = subjects(s);
    if exist(subj.datafile,'file')
       continue; 
       warp=[];
       load(subj.datafile,'warp');
       if isempty(warp)
           warp.to_native = icdm_build_warp_to_native( ...
            subj.icdm_4d, subj.dartel_flow,subj.template);

            warp.to_mni = icdm_build_warp_to_mni( ...
                    subj.icdm_4d, subj.dartel_flow, subj.template);
            save(datafile,'warp','-append','-v7.3');
       end
    else
        fprintf('Data file for subject %s does not exist.', subj.id);

        datafile = fullfile(subj.outpath, sprintf('icdm_subj_%s.mat', subj.id));
        subj.datafile = datafile;

        % Read ICDM 4D
        V  = spm_vol(subj.icdm_4d);
        C4 = spm_read_vols(V);
        thr = getfield_default(opts,'thresh',5);
        sumC = sum(C4,4);
        mask_3d = sumC > thr;
        idx_native = find(mask_3d(:));
        icdm2d = reshape(C4, [], size(C4,4));
        icdm2d = icdm2d(idx_native, :);
        dim_native = V(1).dim;
        V=V(1);
    
        warp.to_native = icdm_build_warp_to_native( ...
                subj.icdm_4d, subj.dartel_flow,subj.template);
    
        warp.to_mni = icdm_build_warp_to_mni( ...
                subj.icdm_4d, subj.dartel_flow, subj.template);
        save(datafile, ...
            'icdm2d','V','idx_native','dim_native','warp','thr','-v7.3');

    end
    subjects(s) = subj;
end
