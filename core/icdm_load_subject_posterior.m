function P = icdm_load_subject_posterior(subj, K)
% ICDM_LOAD_SUBJECT_POSTERIOR  Load index-based VB posterior and rebuild 4-D arrays.
%
%   P = icdm_load_subject_posterior(subj, K)
%
%   Reads the compact index-form posterior ILR map and kappa field saved
%   by icdm_subject_vb, and reconstructs full 4-D volumes suitable for
%   warping or visualisation.
%
%   INPUT
%     subj : subject struct with .outdir
%     K    : number of compositional components
%
%   OUTPUT
%     P : struct with fields
%       .post_ilr : [Nx x Ny x Nz x (K-1)] ILR posterior MAP
%       .mask     : 3-D logical mask
%       .H        : Helmert submatrix
%       .kappa    : [Nx x Ny x Nz] posterior reliability
%
%   Author: Hae-Jeong Park, Ph.D.
%
%   See also ICDM_SUBJECT_VB, ICDM_WARP_TO_MNI

post_file = fullfile(subj.outdir, 'post_ilr_index.mat');
kap_file  = fullfile(subj.outdir, 'kappa_index.mat');

if ~exist(post_file,'file')
    error('[LOAD] Missing posterior: %s', subj.id);
end
if ~exist(kap_file,'file')
    error('[LOAD] Missing kappa: %s', subj.id);
end

S = load(post_file);  % Z0, mask, Nx,Ny,Nz
Kp = K-1;

% reconstruct 4D ILR map
post_ilr = zeros(S.Nx*S.Ny*S.Nz, Kp, 'single');
post_ilr(S.mask(:),:) = S.Z0;
post_ilr = reshape(post_ilr, [S.Nx,S.Ny,S.Nz,Kp]);

P.post_ilr = post_ilr;
P.mask     = S.mask;
P.H        = S.H;

KAP = load(kap_file);
kap = zeros(S.Nx*S.Ny*S.Nz,1,'single');
kap(S.mask(:)) = KAP.kap0;
P.kappa = reshape(kap, [S.Nx,S.Ny,S.Nz]);

end
