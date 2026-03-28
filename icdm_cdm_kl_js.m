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

function [out_kl, out_js, out_nv] = icdm_cdm_kl_js(nii1, nii2, mask_file, out_prefix)
% ICDM_CDM_KL_JS  Compute voxelwise KL and Jensen-Shannon divergence between two CDM NIfTI files.
%
%   [out_kl, out_js, out_nv] = icdm_cdm_kl_js(nii1, nii2, mask_file, out_prefix)
%
%   Loads two ICDM/CDM NIfTI volumes (raw streamline counts or
%   probability maps), normalizes them to the probability simplex,
%   and computes per-voxel KL(p||q) and Jensen-Shannon divergence.
%   Optionally writes a streamline count map (Nv) when inputs are
%   in count mode.  Results are saved as NIfTI files.
%
%   Inputs
%     nii1       : path to first CDM NIfTI (4-D, K volumes).
%     nii2       : path to second CDM NIfTI (same geometry).
%     mask_file  : (optional) path to binary mask NIfTI; auto-computed
%                  from streamline counts if empty.
%     out_prefix : (optional) prefix for output file names (default 'icdm').
%
%   Outputs
%     out_kl : path to the saved KL divergence NIfTI.
%     out_js : path to the saved JS divergence NIfTI.
%     out_nv : path to the saved Nv (streamline count) NIfTI, or '' if
%              inputs were probabilities.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_KL_JS

    %% --------------------------------------------------------
    % defaults
    %% --------------------------------------------------------
    if nargin < 3 || isempty(mask_file)
        mask_file = [];
    end
    if nargin < 4 || isempty(out_prefix)
        out_prefix = 'icdm';
    end

    EPS = 1e-12;
    out_nv = '';

    %% --------------------------------------------------------
    % load NIfTI data
    %% --------------------------------------------------------
    V1 = spm_vol(nii1);
    V2 = spm_vol(nii2);

    % restrict to cortical ROIs
    V1 = V1(1:68);
    V2 = V2(1:68);

    C1 = spm_read_vols(V1);
    C2 = spm_read_vols(V2);

    [X,Y,Z,K] = size(C1);

    %% --------------------------------------------------------
    % detect count vs probability mode
    %% --------------------------------------------------------
    countmode = (max(C1(:)) > 10 || max(C2(:)) > 10);

    %% --------------------------------------------------------
    % mask
    %% --------------------------------------------------------
    if ~isempty(mask_file)
        Vmask = spm_vol(mask_file);
        mask  = spm_read_vols(Vmask) > 0;
    else
        if countmode
            mask = sum(C1,4) > 10 & sum(C2,4) > 10;
        else
            mask = sum(C1,4) > 0.5 & sum(C2,4) > 0.5;
        end
    end

    %% --------------------------------------------------------
    % compute Nv and pi upfront
    %% --------------------------------------------------------
    if countmode
        Nv1 = sum(C1,4);
        Nv2 = sum(C2,4);
        Nvmap = Nv1 + Nv2;

        % normalize to probability
        P1 = C1 ./ max(Nv1, EPS);
        P2 = C2 ./ max(Nv2, EPS);
    else
        P1 = C1;
        P2 = C2;
        Nvmap = [];
    end

    % numerical stabilization (vectorized)
    P1 = max(P1, EPS);
    P2 = max(P2, EPS);

    % renormalize (important after EPS clipping)
    P1 = P1 ./ sum(P1,4);
    P2 = P2 ./ sum(P2,4);

    %% --------------------------------------------------------
    % initialize outputs
    %% --------------------------------------------------------
    KLmap = NaN(X,Y,Z);
    JSmap = NaN(X,Y,Z);

    %% --------------------------------------------------------
    % voxelwise KL / JS (lightweight loop)
    %% --------------------------------------------------------
    for ix = 1:X
        for iy = 1:Y
            for iz = 1:Z

                if ~mask(ix,iy,iz)
                    continue;
                end

                p = squeeze(P1(ix,iy,iz,:));
                q = squeeze(P2(ix,iy,iz,:));

                m = 0.5 * (p + q);

                KLmap(ix,iy,iz) = sum(p .* log(p ./ q));

                JSmap(ix,iy,iz) = ...
                    0.5 * sum(p .* log(p ./ m)) + ...
                    0.5 * sum(q .* log(q ./ m));
            end
        end
    end

    %% --------------------------------------------------------
    % save outputs
    %% --------------------------------------------------------
    Vout = V1(1);
    Vout.dt(1) = 16;  % float32

    if countmode
        out_nv = [out_prefix '_Nv.nii'];
        Vout.fname   = out_nv;
        Vout.descrip = 'Voxelwise streamline count (Nv)';
        spm_write_vol(Vout, Nvmap);
    end

    out_kl = [out_prefix '_KL12.nii'];
    Vout.fname   = out_kl;
    Vout.descrip = 'KL divergence (p || q)';
    spm_write_vol(Vout, KLmap);

    out_js = [out_prefix '_JS.nii'];
    Vout.fname   = out_js;
    Vout.descrip = 'Jensen-Shannon divergence';
    spm_write_vol(Vout, JSmap);

end

