function Ymni = icdm_warp_to_mni(Ynat, warpobj, idx_mni,idx_native,interp)
% ICDM_WARP_TO_MNI  Fast native-to-MNI warping via prebuilt trilinear operator.
%
%   Ymni = icdm_warp_to_mni( Ynative, warpobj, idx_mni, idx_native, interp )
%
%   Fast warping of native-space data to MNI-space data using a
%   prebuilt trilinear warp operator (warpobj.to_mni).
%
%   Used to obtain subject-level results in MNI space:
%       • Posterior ILR mean (y_ilr_native → y_ilr_mni)
%       • Posterior κ-field (kappa_native → kappa_mni)
%       • Any tensor field for group pooling
%
% ----------------------------------------------------------------------
% INPUTS
%
%   Ynative     : Native-space data
%                 Allowed shapes:
%                     [Nnative × D]         (index form)
%                     [Xnat Ynat Znat D]    (full volume)
%
%   warpobj     : Subject’s warp operator containing:
%
%                       warpobj.to_mni.index   [N_mni × 8]
%                       warpobj.to_mni.weight  [N_mni × 8]
%                       warpobj.to_mni.dim_in  = [Xnat Ynat Znat]
%                       warpobj.to_mni.dim_out = [X_mni Y_mni Z_mni]
%
%   idx_mni     : Linear indices of MNI voxels defining the ICDM mask.
%                 Output will be returned as [Nmni × D].
%
%   idx_native  : Linear indices of native WM voxels used in index-form
%                 representation.
%
%   interp      : Interpolation mode:0=nearest neighbor, 1=trilinear (default=1)
%
% ----------------------------------------------------------------------
% OUTPUT
%
%   Ymni        : MNI-space warped data
%                 Always returned as [Nmni × D].
%
% ----------------------------------------------------------------------
% DESCRIPTION
%
%   Steps:
%     1) Convert index-form input into full 4D native volume (if needed).
%     2) Apply fast warp (icdm_apply_warp) using trilinear interpolation.
%     3) Extract only the analysis-space MNI voxels (idx_mni).
%
%   This function **completely replaces** spm_deformations or
%   icdm_warp_4d during subject-VB loops and group-level pooling.
%
%   Accuracy matches SPM DARTEL pull/push within machine precision
%   when used with warp operators built by:
%
%       icdm_build_warp_to_native
%       icdm_build_warp_to_mni
%
% ----------------------------------------------------------------------
% EXAMPLE
%
%   % Warp subject posterior (native → MNI)
%   y_ilr_mni = icdm_warp_to_mni( ...
%       OUT.y_ilr_native, subj.warp, grp.idx_mni, subj.idx_native );
%
%   % Warp κ-field (native → MNI)
%   kappa_mni = icdm_warp_to_mni( ...
%       OUT.kappa_native, subj.warp, grp.idx_mni, subj.idx_native );
%
% ----------------------------------------------------------------------
% NOTE
%
%   For label data, use icdm_warp_to_native with nearest (0).  Native→MNI
%   label warping is usually not required, but can be added if needed.
%
% ----------------------------------------------------------------------
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_WARP_TO_NATIVE, ICDM_APPLY_WARP, ICDM_BUILD_WARP_TO_MNI
% ======================================================================

if nargin <5, interp=1; end
if nargin <4, idx_native = []; end
if nargin <3, idx_mni = []; end
sz=size(Ynat);D=sz(end);
if numel(size(Ynat)) == 2
    P=zeros([warpobj.to_mni.dim_in D],'single');
    for d = 1:D
        p=zeros(warpobj.to_mni.dim_in); 
        p(idx_native)=Ynat(:,d);
        P(:,:,:,d) = p;
    end
    Ynat = P;
end

nnflag=0;
if interp==0 && D==1 % nearest neighbor
    Kmax = max(Ynat(:));
    dim = size(Ynat);
    P = zeros([dim Kmax],'single');
    for k = 1:Kmax
        P(:,:,:,k) = single(Ynat == k);
    end
    Ynat=P;
    nnflag=1;
end

Ymni = icdm_apply_warp(Ynat,warpobj.to_mni);

if nnflag, [~, Ymni] = max(Ymni, [], 4); end

if nargin>=3 && ~isempty(idx_mni)
    D=size(Ymni,4);
    Ymni= reshape(Ymni,[],D);
    Ymni = Ymni(idx_mni,:);
end
end
