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

function icdm_fs_mgz_to_vep(mgzpath,veppath)
% ICDM_FS_MGZ_TO_VEP  Convert FreeSurfer MGZ parcellations to VEP NIfTI.
%
%   icdm_fs_mgz_to_vep(mgzpath, veppath)
%
%   Converts standard FreeSurfer MGZ files (aparc+aseg, orig, brain) from
%   the subject's mri/ directory into RAS-oriented NIfTI files, and
%   generates the VEP-166 parcellation from aparc+aseg.vep.mgz using
%   mnet_generate_vep166_from_aparc.
%
%   INPUT
%     mgzpath : path to FreeSurfer subject directory (containing mri/)
%     veppath : output directory for NIfTI files
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_COMPOSE_SUBJECT
    mgz_names={'aparc+aseg.mgz',...
            'aparc+aseg.vep.mgz',...
            'orig.mgz','brain.mgz'};
    for i=1:length(mgz_names)
        mgzfile=fullfile(mgzpath,'mri',mgz_names{i});
        niifile=fullfile(veppath,strrep(mgz_names{i},'.mgz','.nii'));
        if exist(mgzfile,'file') && ~exist(niifile,'file')
            cmd=sprintf('!mri_convert --out_orientation RAS %s %s',mgzfile,niifile);
            eval(cmd);
            fprintf('Converted %s to %s\n', mgzfile, niifile);  
            % where %s are the mgzfile and niifile paths
        end
        if strcmp(mgz_names{i},'aparc+aseg.vep.mgz')
            niifile_vep=fullfile(veppath,'VEP_166.nii');
            if exist(niifile,'file') && ~exist(niifile_vep,'file')
                mnet_generate_vep166_from_aparc(niifile, niifile_vep);
                fprintf('Converted %s to %s\n', niifile, niifile_vep);  
            end
        end
    end
end

