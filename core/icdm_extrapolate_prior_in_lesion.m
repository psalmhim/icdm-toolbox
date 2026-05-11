function [MU_out, KAP_out, affected_mask] = icdm_extrapolate_prior_in_lesion( ...
    MU_nat, KAP_nat, lesion_vol, idx_native, dim_native, varargin)
% ICDM_EXTRAPOLATE_PRIOR_IN_LESION  Displacement-aware prior mapping for
%   lesion/tumor-affected brains.
%
%   When a tumor distorts white-matter anatomy through mass effect, the
%   standard deformation-field-based group prior mapping becomes unreliable
%   in and around the lesion. This function replaces the group prior at
%   affected voxels with values extrapolated from nearby healthy tissue
%   using Gaussian-weighted spatial inpainting.
%
%   Conceptually, this provides a displacement-corrected prior: instead of
%   using the (distorted) MNI coordinate for a tumor-adjacent voxel, the
%   prior is interpolated from healthy neighbors whose MNI coordinates are
%   reliable, approximating the prior that would have been obtained if the
%   tissue were in its pre-displacement position.
%
%   [MU_out, KAP_out, affected] = icdm_extrapolate_prior_in_lesion( ...
%       MU_nat, KAP_nat, lesion_vol, idx_native, dim_native, ...)
%
%   INPUTS
%     MU_nat      [Nnat x K1]  Group prior mean in native space (from warp)
%     KAP_nat     [Nnat x 1]   Group prior reliability in native space
%     lesion_vol  [X x Y x Z]  Lesion probability map in native space
%                               (values in [0,1]; 0 = healthy, 1 = lesion)
%     idx_native  [Nnat x 1]   Linear indices of WM mask voxels in native vol
%     dim_native  [1 x 3]      Native volume dimensions [X Y Z]
%
%   OPTIONAL NAME-VALUE PAIRS
%     'LesionThreshold'   scalar   Threshold on lesion_vol to define lesion
%                                  core (default: 0.3)
%     'DilationRadius'    scalar   Dilation radius in voxels to define the
%                                  mass-effect affected zone around the
%                                  lesion (default: 5)
%     'SigmaVox'          scalar   Gaussian kernel sigma in voxels for
%                                  spatial extrapolation (default: 4)
%     'KappaFloor'        scalar   Minimum kappa after correction
%                                  (default: 0.1)
%     'AttenuateKappa'    logical  Whether to attenuate kappa in the
%                                  affected zone (default: true)
%     'KappaAttenuation'  scalar   Factor by which kappa is reduced in the
%                                  affected zone (default: 0.5)
%
%   OUTPUTS
%     MU_out      [Nnat x K1]  Corrected prior mean
%     KAP_out     [Nnat x 1]   Corrected prior reliability
%     affected    [X x Y x Z]  Binary mask of the affected zone
%
%   ALGORITHM
%     1. Threshold lesion_vol to create a binary lesion core mask
%     2. Dilate the core by DilationRadius voxels to define the
%        mass-effect affected zone (where the deformation field is
%        unreliable)
%     3. Define healthy voxels = WM mask AND NOT affected zone
%     4. Embed MU and KAP into 3D volumes, masked to healthy voxels only
%     5. Apply Gaussian smoothing (sigma = SigmaVox) to both the
%        healthy-masked data and the healthy mask itself
%     6. Normalise: corrected = smoothed_data / smoothed_mask
%        This yields a distance-weighted average from healthy voxels
%     7. Replace values at affected WM voxels with the extrapolated values
%     8. Optionally attenuate kappa in the affected zone to reflect
%        reduced confidence in the extrapolated prior
%
%   REFERENCE
%     Park et al., "Reliability-Aware Hierarchical Bayesian Inference for
%     Voxelwise Compositional Connectivity Mapping in Diffusion MRI"
%
%   See also ICDM_SUBJECT_VB, ICDM_WARP_TO_NATIVE

% ---- parse options ----
p = inputParser;
addParameter(p, 'LesionThreshold', 0.3, @isscalar);
addParameter(p, 'DilationRadius', 5, @isscalar);
addParameter(p, 'SigmaVox', 4, @isscalar);
addParameter(p, 'KappaFloor', 0.1, @isscalar);
addParameter(p, 'AttenuateKappa', true, @islogical);
addParameter(p, 'KappaAttenuation', 0.5, @isscalar);
parse(p, varargin{:});
o = p.Results;

[Nnat, K1] = size(MU_nat);
Xn = dim_native(1); Yn = dim_native(2); Zn = dim_native(3);

% ---- 1. Binary lesion core ----
lesion_core = lesion_vol > o.LesionThreshold;

% ---- 2. Dilate to define affected zone ----
if o.DilationRadius > 0
    se = strel('sphere', round(o.DilationRadius));
    affected_mask = imdilate(lesion_core, se);
else
    affected_mask = lesion_core;
end

% ---- 3. Build WM mask volume and healthy mask ----
wm_mask = false(Xn, Yn, Zn);
wm_mask(idx_native) = true;

healthy_mask = wm_mask & ~affected_mask;
affected_wm  = wm_mask &  affected_mask;

% Check how many voxels are affected
n_affected = nnz(affected_wm);
if n_affected == 0
    MU_out  = MU_nat;
    KAP_out = KAP_nat;
    fprintf('  [Lesion prior] No WM voxels in affected zone; no correction applied.\n');
    return;
end
fprintf('  [Lesion prior] %d / %d WM voxels in affected zone (%.1f%%)\n', ...
    n_affected, Nnat, 100*n_affected/Nnat);

% ---- 4. Embed into 3D volumes (healthy voxels only) ----
healthy_vol = single(healthy_mask);

MU_vol  = zeros(Xn, Yn, Zn, K1, 'single');
KAP_vol = zeros(Xn, Yn, Zn, 'single');

% Map masked indices back to volume
for d = 1:K1
    tmp = zeros(Xn, Yn, Zn, 'single');
    tmp(idx_native) = MU_nat(:, d);
    tmp(~healthy_mask) = 0;              % zero out non-healthy
    MU_vol(:,:,:,d) = tmp;
end

tmp = zeros(Xn, Yn, Zn, 'single');
tmp(idx_native) = KAP_nat;
tmp(~healthy_mask) = 0;
KAP_vol = tmp;

% ---- 5. Gaussian smoothing ----
sigma = o.SigmaVox;
ksz   = ceil(3 * sigma) * 2 + 1;       % kernel support: 3σ each side

% Smooth the healthy mask (normalisation denominator)
healthy_smooth = imgaussfilt3(healthy_vol, sigma, 'FilterSize', ksz);
healthy_smooth = max(healthy_smooth, eps('single'));

% Smooth MU (each ILR dimension separately)
MU_smooth = zeros(size(MU_vol), 'single');
for d = 1:K1
    MU_smooth(:,:,:,d) = imgaussfilt3(MU_vol(:,:,:,d), sigma, 'FilterSize', ksz);
end

% Smooth KAP
KAP_smooth = imgaussfilt3(KAP_vol, sigma, 'FilterSize', ksz);

% ---- 6. Normalise: distance-weighted average from healthy voxels ----
for d = 1:K1
    MU_smooth(:,:,:,d) = MU_smooth(:,:,:,d) ./ healthy_smooth;
end
KAP_smooth = KAP_smooth ./ healthy_smooth;

% ---- 7. Replace values at affected WM voxels ----
MU_out  = MU_nat;
KAP_out = KAP_nat;

% Build lookup: which masked-index voxels are in the affected zone?
affected_lin = find(affected_wm(:));
[~, ia, ib] = intersect(idx_native, affected_lin);
% ia = indices into MU_nat / KAP_nat
% ib = indices into affected_lin (not needed further)

for d = 1:K1
    slice = MU_smooth(:,:,:,d);
    MU_out(ia, d) = slice(affected_lin(ib));
end
KAP_out(ia) = KAP_smooth(affected_lin(ib));

% ---- 8. Attenuate kappa in affected zone ----
if o.AttenuateKappa
    % Graded attenuation: stronger near lesion core, weaker at periphery
    dist_to_core = bwdist(lesion_core);
    attenuation = zeros(Xn, Yn, Zn, 'single');
    attenuation(affected_mask) = o.KappaAttenuation + ...
        (1 - o.KappaAttenuation) * min(dist_to_core(affected_mask) / o.DilationRadius, 1);
    % attenuation = KappaAttenuation at core boundary, → 1.0 at edge of affected zone
    KAP_out(ia) = KAP_out(ia) .* attenuation(affected_lin(ib));
end

% Floor
KAP_out = max(KAP_out, o.KappaFloor);

fprintf('  [Lesion prior] Extrapolation complete. sigma=%.1f, dilation=%d vox\n', ...
    sigma, round(o.DilationRadius));

end
