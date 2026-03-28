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

function Ynat = icdm_warp_to_native(Ymni, warpobj, idx_native, idx_mni,interp)
% ICDM_WARP_TO_NATIVE  Fast MNI-to-native warping via prebuilt trilinear operator.
%
%   Ynat = icdm_warp_to_native( Ymni, warpobj, idx_native, idx_mni, interp )
%
%   Fast warping of MNI-space data to native space using a prebuilt
%   trilinear warp operator (warpobj.to_native).
%
%   Supports:
%     - index-form Ymni  [Nmni x D]  reshaped to full MNI grid
%     - volume Ymni      [X Y Z D]
%     - label warping    (nearest-neighbor via trilinear argmax)
%     - masked return with idx_native
%
%   Inputs
%     Ymni       : MNI-space data, [Nmni x D] or [X Y Z D].
%     warpobj    : warp struct containing warpobj.to_native.
%     idx_native : (optional) linear indices for masked native output.
%     idx_mni    : (optional) linear indices for index-form MNI input.
%     interp     : 0 = nearest neighbor, 1 = trilinear (default 1).
%
%   Output
%     Ynat : native-space warped data, [Nnative x D] if idx_native given.
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_WARP_TO_MNI, ICDM_APPLY_WARP, ICDM_BUILD_WARP_TO_NATIVE
if nargin <5, interp=1; end
if nargin <4, idx_mni = []; end
if nargin <3, idx_native = []; end

sz=size(Ymni);D=sz(end);
if numel(size(Ymni)) == 2
    P=zeros([warpobj.to_native.dim_in D],'single');
    for d = 1:D
        p=zeros(warpobj.to_native.dim_in); p(idx_mni)=Ymni(:,d);
        P(:,:,:,d) = p;
    end
    Ymni = P;
end

nnflag=0;
if interp==0 && D==1 % nearest neighbor
    Kmax = max(Ymni(:));
    dim = size(Ymni);
    P = zeros([dim Kmax],'single');
    for k = 1:Kmax
        P(:,:,:,k) = single(Ymni == k);
    end
    Ymni=P;
    nnflag=1;
end
Ynat = icdm_apply_warp(Ymni,warpobj.to_native);


if nnflag, [~, Ynat] = max(Ynat, [], 4); end

if nargin>=3
    D=size(Ynat,4);
    Ynat= reshape(Ynat,[],D);
    Ynat = Ynat(idx_native,:);
end
end
