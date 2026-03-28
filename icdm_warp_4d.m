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

function Yout = icdm_warp_4d(Yin, nativefile, def_field, templatefile, direction, interp, temp_dir)
% ICDM_WARP_4D  Warp a 4-D array between native and MNI space via SPM12.
%
%   Yout = icdm_warp_4d(Yin, nativefile, def_field, templatefile, direction, interp, temp_dir)
%
%   Frame-by-frame warping of a 4-D volume through the DARTEL deformation
%   field using spm_deformations.  Each 3-D frame is written to a
%   temporary NIfTI, warped, read back, and cleaned up.
%
%   INPUT
%     Yin          : [X x Y x Z x D] 4-D input array
%     nativefile   : native-space reference (path or spm_vol struct)
%     def_field    : path to DARTEL deformation field (y_*.nii)
%     templatefile : path to DARTEL template NIfTI
%     direction    : +1 = native-to-MNI, -1 = MNI-to-native
%     interp       : interpolation order (default 1 = trilinear)
%     temp_dir     : directory for temporary files (default pwd)
%
%   OUTPUT
%     Yout : warped 4-D single array in the target space
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_APPLY_WARP, GENERIC_WARP_4D, SPM_DEFORMATIONS

warning('off','nifti:hdr');
warning('off','nifti:datatype');
warning('off','SPM:hdr:tooBig');

if nargin<6 || isempty(interp), interp = 1; end
if nargin<7, temp_dir = ''; end
% ----------------------------------------------------------
% Temp directory
% ----------------------------------------------------------
persistent TEMP_DIR
if isempty(TEMP_DIR)
    if exist('temp_dir','var') && ~isempty(temp_dir)
        TEMP_DIR = temp_dir;
    else
        TEMP_DIR = pwd;
    end
end

if isempty(def_field) || exist(def_field,'file')~=2
    warning('[DARTEL] deformation missing -> pass-through');
    Yout = Yin;
    return;
end

if direction>0 && (isempty(templatefile) || exist(templatefile,'file')~=2)
    error('icdm_warp_4d: templatefile required for forward direction.');
end

Yin = single(Yin);
[~,~,~,D] = size(Yin);

% ----------------------------------------------------------
% Helper: write temporary volume
% ----------------------------------------------------------
function Vi = write_tmp_vol(vol3d, Vlike, tag)
    uuid = char(java.util.UUID.randomUUID);
    baseDir = TEMP_DIR;
    fname = fullfile(baseDir, ['warp_' tag '_' uuid '.nii']);
    Vi = Vlike;
    Vi.fname = fname;
    Vi.dt = [spm_type('float32') 0];
    Vi.pinfo = [1;0;0];
    Vi.n = [1 1];
    Vi.descrip = 'icdm_warp_4d tmp';
    sz = size(vol3d);
    if numel(sz) < 3, sz = [sz, ones(1,3-numel(sz))]; end
    Vi.dim(1:3) = sz(1:3);
    spm_write_vol(Vi, single(vol3d));
end

% ----------------------------------------------------------
% Define reference for SPM job
% ----------------------------------------------------------
if direction > 0
    % native → MNI
    if ischar(nativefile)
        vnative = spm_vol(nativefile);
    else
        vnative = nativefile;
    end
    Vin_like = vnative(1);
else
    % MNI → native
    if ischar(nativefile)
        vnative = spm_vol(nativefile);
    else
        vnative = nativefile;
    end
    Vtpl = spm_vol(templatefile);
    Vin_like = Vtpl(1);
end

% ----------------------------------------------------------
% Loop over each frame
% ----------------------------------------------------------
Yout = [];
warped_files = strings(D,1);

try
    for d = 1:D
        
        Vi = write_tmp_vol(Yin(:,:,:,d), Vin_like, sprintf('d%03d',d));

        job = struct();
        job.comp = {};

        if direction > 0
            % native → MNI
            job.comp{1}.def = {def_field};
            job.space       = {templatefile};
        else
            % MNI → native
            job.comp{1}.inv.comp{1}.def = {def_field};
            job.comp{1}.inv.space       = {vnative(1).fname};
        end

        job.out{1}.pull.fnames          = {Vi.fname};
        job.out{1}.pull.savedir.saveusr = {TEMP_DIR};
        job.out{1}.pull.interp          = interp;
        job.out{1}.pull.mask            = 0;
        job.out{1}.pull.fwhm            = [0 0 0];

        out = spm_deformations(job);
        if isempty(out.warped)
            error('spm_deformations returned empty output.');
        end

        fout = out.warped{1};
        warped_files(d) = string(fout);

        Yo = spm_read_vols(spm_vol(fout));

        if isempty(Yout)
            [Xo,Yo_,Zo] = size(Yo);
            Yout = zeros(Xo,Yo_,Zo,D,'single');
        end

        Yout(:,:,:,d) = single(Yo);

        if exist(Vi.fname,'file'), try, delete(Vi.fname); end, end
        if exist(fout,'file'), try, delete(fout); end, end
    end

catch ME

    % cleanup
    tmpfiles = dir(fullfile(TEMP_DIR, 'warp_*'));
    for tt = 1:numel(tmpfiles)
        try, delete(fullfile(TEMP_DIR, tmpfiles(tt).name)); end
    end

    rethrow(ME);
end

% Final cleanup
for dd = 1:numel(warped_files)
    f = warped_files(dd);
    if ~isempty(f) && exist(f,'file')
        try, delete(f); end
    end
end

Yout = single(Yout);

end