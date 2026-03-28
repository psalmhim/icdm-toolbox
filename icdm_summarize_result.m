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

% ICDM_SUMMARIZE_RESULT  Interactive script to visualize group iCDM iteration outputs.
%
%   Loads NIfTI volumes (kappa, mu_ilr, pi, tau2, beta) produced by
%   icdm_result_to_nii and displays coronal tile mosaics and multi-slice
%   overlays for each ILR dimension, pi parcel, and beta field.
%   Requires a T1 reference image (rbMNI152_mri.nii) in the output
%   directory.
%
%   See also ICDM_RESULT_TO_NII, ICDM_DISPLAY_ALL
% Author: Hae-Jeong Park, Ph.D.

%function icdm_summarize_result(outdir,itername)

%% ===== USER CONFIGURATION (edit before running) =====
% outdir   = '/path/to/your/ICDM_output';
% itername = '003';
error('Edit the USER CONFIGURATION section above before running this script.');

outdir='/home/hjpark/data/ICDM_desikan68';
itername='003';
if ~exist('KAPPA_mni','var')
    fn = fullfile(outdir,[itername '_group_kappa.nii']);
    v=spm_vol(fn);
    KAPPA_mni = spm_read_vols(v);
    
    fn = fullfile(outdir,[itername '_group_mu_ilr.nii']);
    v=spm_vol(fn);
    MU_ilr_mni = spm_read_vols(v);
    
    fn = fullfile(outdir,[itername '_group_pi.nii']);
    v=spm_vol(fn);
    PI_mni = spm_read_vols(v);
    
    fn = fullfile(outdir,[itername '_group_tau2_ilr.nii']);
    v=spm_vol(fn);
    TAU_mni = spm_read_vols(v);
    
    fn = fullfile(outdir,[itername '_group_beta1.nii']);
    v=spm_vol(fn);
    Beta1_mni = spm_read_vols(v);
    
    fn = fullfile(outdir,[itername '_group_beta2.nii']);
    v=spm_vol(fn);
    Beta2_mni = spm_read_vols(v);
    
    vt=spm_vol(fullfile(outdir,'rbMNI152_mri.nii'));
    T=spm_read_vols(vt);
end

figure;
show_coronal_tiles(MU_ilr_mni,T);

figure;
show_multi_slices(MU_ilr_mni,T,1);
figure;
show_coronal_tiles(PI_mni,T);
figure;
show_multi_slices(PI_mni,T, 1);




%% ---------------------------------------------------------
% Montage view
%% ---------------------------------------------------------
% figure('Name','Montage of Kappa','Color','w');
% montage(PI_mni,'DisplayRange',[]);
% title('Group Reliability (Kappa)');
% 
% fprintf('Visualization complete.\n');

function show_coronal_tiles(R,T,ratio)
if nargin<3, ratio=0.0; end
K1=size(R,4);
dim=size(R);hdim=round(dim/2);
Rs=squeeze(R(:,hdim(2),:,:));
Ts=squeeze(T(:,hdim(2),:));
mR=max(abs(Rs(:)))*0.9; 
thr=mR*ratio;
t = tiledlayout(ceil(sqrt(K1)), ceil(sqrt(K1)));
t.TileSpacing = 'none';
t.Padding = 'none';
for d = 1:K1
    nexttile;
    rR = corient(Rs(:,:,d));
    rT = corient(Ts); 
    monet_overlay(rT,gray(256),[],rR,[-thr thr],jet(256),[-mR mR],0.2);
    axis image off;
end
colormap coolwarm;
end

function show_multi_slices(R,T,region,cmap,ratio)
if nargin<4,cmap=[]; end
if nargin<5, ratio=0.1; end
R = squeeze(R(:,:,:,region));
[m,mid]=max(R(:)); [x,y,z]=ind2sub(size(R),mid);
mR=max(abs(R(:)))*0.9; 
thr=mR*ratio;
if isempty(cmap)
    if min(R(:))<0
     cmap=coolwarm(256);
    else
        cmap=hot(256);
    end
end
subplot(1,3,1);
rR = corient(R(:,:,z));
rT = corient(T(:,:,z));
monet_overlay(rT,gray(256),[],rR,[-thr thr],cmap,[-mR mR],0.2);
axis image off;
title(sprintf('Region %d: Axial', region));

subplot(1,3,2);
rR = corient(R(:,y,:));
rT = corient(T(:,y,:));
monet_overlay(rT,gray(256),[],rR,[-thr thr],cmap,[-mR mR],0.2);
axis image off;
title('coronal');

subplot(1,3,3);
rR = corient(R(x,:,:));
rT = corient(T(x,:,:));
monet_overlay(rT,gray(256),[],rR,[-thr thr],cmap,[-mR mR],0.2);
axis image off;
title('Sagittal');

end


function img=corient(img)
    img = rot90(squeeze(img), -1);
end