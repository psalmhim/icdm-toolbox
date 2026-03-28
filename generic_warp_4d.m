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

function Yout = generic_warp_4d(Yin, Vref_native, def_field, template, direction, interp)
% GENERIC_WARP_4D  Warp a 4-D array between native and MNI space via SPM12.
%
%   Yout = generic_warp_4d(Yin, Vref_native, def_field, template, direction, interp)
%
%   Wrapper around spm_deformations that writes each 3-D frame of Yin to a
%   temporary NIfTI, applies the DARTEL deformation field, reads the warped
%   result, and stacks the frames back into a 4-D output.  Temporary files
%   are cleaned up automatically even on error.
%
%   Inputs
%     Yin         : 4-D single array [X x Y x Z x D].
%     Vref_native : spm_vol header for the native-space reference.
%     def_field   : path to DARTEL deformation field (y_*.nii).
%     template    : path to DARTEL template NIfTI.
%     direction   : +1 = native-to-MNI (forward pull),
%                   -1 = MNI-to-native (inverse pull).
%     interp      : interpolation order (default 1 = trilinear).
%
%   Output
%     Yout : warped 4-D single array in the target space.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_WARP_4D, SPM_DEFORMATIONS

    if nargin<6 || isempty(interp), interp = 1; end
    if isempty(def_field) || exist(def_field,'file')~=2
        warning('[DARTEL] deformation missing -> pass-through'); Yout = Yin; return;
    end
    if direction>0 && (isempty(template) || exist(template,'file')~=2)
        error('generic_warp_4d: template required for forward pulls.');
    end
    Yin = single(Yin);
    [~,~,~,D] = size(Yin);

    function Vi = write_tmp_vol(vol3d, Vlike, tag)
        if nargin<3, tag=''; end
        baseDir = fileparts(char(Vlike.fname));
        Vi = Vlike; Vi.fname = fullfile(baseDir, ['tmp_' tag '_' char(java.util.UUID.randomUUID) '.nii']);
        [pth,~] = fileparts(Vi.fname); if ~exist(pth,'dir'), mkdir(pth); end
        Vi.dt=[spm_type('float32') 0]; Vi.pinfo=[1;0;0]; Vi.n=[1 1]; Vi.descrip='generic_warp_4d tmp';
        sz = size(vol3d); if numel(sz)<3, sz = [sz, ones(1,3-numel(sz))]; end
        Vi.dim(1:3) = sz(1:3);
        spm_write_vol(Vi, single(vol3d));
    end

    if direction > 0
        Vin_like      = Vref_native(1);
        target_spaceF = char(template);
    else
        Vtpl = spm_vol(template);
        Vin_like      = Vtpl(1);
        target_spaceF = char(Vref_native(1).fname);
    end

    Yout = []; warped_files = strings(D,1);
    try
        for d = 1:D
            Vi = write_tmp_vol(Yin(:,:,:,d), Vin_like, sprintf('d%03d',d));
            job = struct(); job.comp = {};
            if direction > 0
                job.comp{1}.def = {def_field};
                job.space       = {target_spaceF};
            else
                job.comp{1}.inv.comp{1}.def = {def_field};
                job.comp{1}.inv.space       = {target_spaceF};
            end
            job.out{1}.pull.fnames          = {Vi.fname};
            job.out{1}.pull.savedir.saveusr = {fileparts(Vi.fname)};
            job.out{1}.pull.interp          = interp;
            job.out{1}.pull.mask            = 0;
            job.out{1}.pull.fwhm            = [0 0 0];
            out = spm_deformations(job);
            if isempty(out.warped), error('spm_deformations returned no warped output.'); end
            fout = out.warped{1}; warped_files(d) = string(fout);
            Yo = spm_read_vols(spm_vol(fout));
            if isempty(Yout), [Xo,Yo_,Zo] = size(Yo); Yout = zeros(Xo,Yo_,Zo,D,'single'); end
            Yout(:,:,:,d) = single(Yo);
            try, delete(Vi.fname); catch, end
        end
    catch ME
        for dd=1:numel(warped_files)
            f = char(warped_files(dd));
            if ~isempty(f) && exist(f,'file'), try, delete(f); catch, end, end
        end
        rethrow(ME);
    end
    for dd=1:numel(warped_files)
        f = char(warped_files(dd));
        if ~isempty(f) && exist(f,'file'), try, delete(f); catch, end, end
    end
end
