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

function [MU2, K2] = spatial_smooth_prior(MU, K, mask3d, dimN, sopts)
% SPATIAL_SMOOTH_PRIOR  Spatial smoothing of group-level ILR prior.
%
%   [MU2, K2] = spatial_smooth_prior(MU, K, mask3d, dimN, sopts)
%
%   Applies one-pass neighbourhood smoothing to the group prior mean and
%   precision in ILR space (manuscript Eqs. spatial_smoothing,
%   spatial_smoothing_precision). Each neighbour contributes weight eta
%   to the effective sample size:
%
%     mu_v <- (kappa_v * mu_v + eta * sum_{u in N(v)} mu_u) /
%             (kappa_v + eta * |N(v)|)
%     kappa_v <- kappa_v + eta * |N(v)|
%
%   INPUT
%     MU     : [Nv x (K-1)] group prior mean in ILR space
%     K      : [Nv x 1] group prior precision (kappa)
%     mask3d : [X x Y x Z] logical brain mask
%     dimN   : [X Y Z] volume dimensions
%     sopts  : struct with .lambda (smoothing weight eta) and
%              .neighborhood (6 or 26 connectivity)
%
%   OUTPUT
%     MU2 : [Nv x (K-1)] spatially smoothed prior mean
%     K2  : [Nv x 1] updated precision (increased by neighbour count)
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_UPDATE_GROUP_PRIOR

lambda = sopts.lambda;
if lambda <= 0
    MU2 = MU;
    K2 = K;
    return;
end

[Xn,Yn,Zn] = deal(dimN(1),dimN(2),dimN(3));
Nm = size(MU,2);

% create linear index → neighbor list
if sopts.neighborhood == 26
    nbr = get_neighbors_26(mask3d, Xn,Yn,Zn);
else
    nbr = get_neighbors_6(mask3d, Xn,Yn,Zn);
end

numV = numel(K);
MU2 = MU;
K2  = K;

for v = 1:numV
    nv = nbr{v};
    if isempty(nv), continue; end
    deg = numel(nv);

    MU_neighbors = MU(nv,:);     % deg × K1

    num = K(v)*MU(v,:) + lambda*sum(MU_neighbors,1);
    den = K(v) + lambda*deg;

    MU2(v,:) = num./den;
    K2(v)    = den;
end

end


%% ========================================================================
%% NEIGHBORHOOD FUNCTIONS
%% ========================================================================
function nbr = get_neighbors_6(mask, X,Y,Z)
idx = find(mask(:));
nbr = cell(numel(idx),1);

% shift offsets for 6-neighborhood
offsets = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];

LUT = zeros(numel(mask),1,'uint32');
LUT(idx) = 1:numel(idx);

for i = 1:numel(idx)
    lin = idx(i);
    [x,y,z] = ind2sub([X Y Z], lin);

    neigh_lin = [];

    for k=1:6
        xn = x+offsets(k,1);
        yn = y+offsets(k,2);
        zn = z+offsets(k,3);
        if xn>=1 && xn<=X && yn>=1 && yn<=Y && zn>=1 && zn<=Z
            j = sub2ind([X Y Z], xn,yn,zn);
            if mask(j)
                neigh_lin(end+1) = LUT(j);
            end
        end
    end
    nbr{i} = neigh_lin;
end

end


function nbr = get_neighbors_26(mask, X,Y,Z)
idx = find(mask(:));
nbr = cell(numel(idx),1);

LUT = zeros(numel(mask),1,'uint32');
LUT(idx) = 1:numel(idx);

for i=1:numel(idx)
    lin=idx(i); [x,y,z]=ind2sub([X Y Z],lin);
    neigh=[];
    for dx=-1:1
        for dy=-1:1
            for dz=-1:1
                if dx==0 && dy==0 && dz==0
                    continue;
                end
                xn=x+dx; yn=y+dy; zn=z+dz;
                if xn>=1 && xn<=X && yn>=1 && yn<=Y && zn>=1 && zn<=Z
                    j=sub2ind([X Y Z],xn,yn,zn);
                    if mask(j)
                        neigh(end+1)=LUT(j);
                    end
                end
            end
        end
    end
    nbr{i}=neigh;
end

end

