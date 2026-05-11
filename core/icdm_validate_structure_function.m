function val = icdm_validate_structure_function(pi_struct, fcdm, subj, opts)
% ICDM_VALIDATE_STRUCTURE_FUNCTION  Compare structural iCDM with functional CDM.
%
%   val = icdm_validate_structure_function(pi_struct, fcdm, subj, opts)
%
%   Evaluates structure-function agreement by computing voxelwise
%   similarity between the structural connectivity distribution (iCDM)
%   and the functional connectivity distribution (fCDM from rs-fMRI).
%
%   For tumor validation, results are stratified by distance from the
%   lesion to test whether displacement-aware prior correction improves
%   structure-function agreement in peritumoral tissue.
%
%   INPUT
%     pi_struct : [Nnat x K] structural connectivity distribution
%                 (from iCDM: ilr_inverse(y_ilr, H))
%     fcdm      : struct from icdm_compute_functional_cdm
%                   .pi_func    [Nnat x K]
%                   .idx_native [Nnat x 1]
%                   .dim_native [1 x 3]
%     subj      : subject struct
%                   .Lesion     [X x Y x Z] lesion mask (optional)
%     opts      : options struct
%       .validate.dist_bins   : distance bin edges in voxels
%                               (default [0 5 10 20 Inf])
%       .validate.exclude_lesion : exclude lesion core from metrics
%                               (default true)
%
%   OUTPUT
%     val : struct with fields
%       .js_voxelwise   : [Nnat x 1] Jensen-Shannon divergence per voxel
%       .cosine_voxelwise : [Nnat x 1] cosine similarity per voxel
%       .spearman_voxelwise : [Nnat x 1] Spearman rho per voxel
%       .dist_to_lesion : [Nnat x 1] distance to lesion boundary (voxels)
%                         (NaN if no lesion)
%       .bin_labels     : cell array of bin labels
%       .js_by_dist     : [Nbins x 1] median JS per distance bin
%       .cosine_by_dist : [Nbins x 1] median cosine per distance bin
%       .spearman_by_dist : [Nbins x 1] median Spearman per distance bin
%       .n_by_dist      : [Nbins x 1] voxel count per bin
%       .summary        : struct with overall median JS, cosine, spearman
%
%   See also ICDM_COMPUTE_FUNCTIONAL_CDM, ICDM_SUBJECT_VB

fprintf('[Validate] Structure-function comparison\n');

% ---- defaults ----
val_opts = getfield_default(opts, 'validate', struct());
dist_bins = getfield_default(val_opts, 'dist_bins', [0 5 10 20 Inf]);
exclude_lesion = getfield_default(val_opts, 'exclude_lesion', true);

pi_func = fcdm.pi_func;
idx_native = fcdm.idx_native;
dim_native = fcdm.dim_native;
[Nnat, K] = size(pi_struct);

assert(size(pi_func,1) == Nnat, 'Structural and functional CDM size mismatch');
assert(size(pi_func,2) == K,   'Structural and functional CDM dimension mismatch');

% ---- 1. Compute voxelwise metrics ----
fprintf('  Computing voxelwise metrics (%d voxels)...\n', Nnat);

js_vox   = zeros(Nnat, 1, 'single');
cos_vox  = zeros(Nnat, 1, 'single');
rho_vox  = zeros(Nnat, 1, 'single');

% Floor to avoid log(0)
eps_floor = 1e-10;
ps = max(double(pi_struct), eps_floor);
pf = max(double(pi_func),   eps_floor);
% Re-normalise after floor
ps = ps ./ sum(ps, 2);
pf = pf ./ sum(pf, 2);

parfor v = 1:Nnat
    p = ps(v,:);
    q = pf(v,:);

    % Jensen-Shannon divergence
    m = 0.5 * (p + q);
    js_vox(v) = 0.5 * sum(p .* log(p ./ m)) + 0.5 * sum(q .* log(q ./ m));

    % Cosine similarity
    cos_vox(v) = sum(p .* q) / (norm(p) * norm(q) + 1e-12);

    % Spearman rank correlation
    rho_vox(v) = corr(p', q', 'Type', 'Spearman');
end

js_vox(~isfinite(js_vox))   = NaN;
cos_vox(~isfinite(cos_vox)) = NaN;
rho_vox(~isfinite(rho_vox)) = NaN;

% ---- 2. Distance to lesion ----
has_lesion = isfield(subj, 'Lesion') && ~isempty(subj.Lesion);
dist_to_lesion = NaN(Nnat, 1, 'single');

if has_lesion
    [Xn, Yn, Zn] = deal(dim_native(1), dim_native(2), dim_native(3));
    lesion_core = subj.Lesion > 0.5;
    dist_vol = single(bwdist(lesion_core));
    dist_to_lesion = dist_vol(idx_native);

    % Optionally exclude lesion core voxels
    if exclude_lesion
        in_lesion = dist_to_lesion < 0.5;  % inside core
        js_vox(in_lesion)  = NaN;
        cos_vox(in_lesion) = NaN;
        rho_vox(in_lesion) = NaN;
        fprintf('  Excluded %d lesion-core voxels from metrics\n', sum(in_lesion));
    end
end

% ---- 3. Stratify by distance bins ----
Nbins = numel(dist_bins) - 1;
bin_labels  = cell(Nbins, 1);
js_by_dist  = NaN(Nbins, 1);
cos_by_dist = NaN(Nbins, 1);
rho_by_dist = NaN(Nbins, 1);
n_by_dist   = zeros(Nbins, 1);

if has_lesion
    fprintf('  Distance bins (voxels): ');
    for b = 1:Nbins
        lo = dist_bins(b);
        hi = dist_bins(b+1);
        in_bin = dist_to_lesion >= lo & dist_to_lesion < hi;
        n_by_dist(b)   = sum(in_bin & isfinite(js_vox));
        js_by_dist(b)  = nanmedian(js_vox(in_bin));
        cos_by_dist(b) = nanmedian(cos_vox(in_bin));
        rho_by_dist(b) = nanmedian(rho_vox(in_bin));

        if isinf(hi)
            bin_labels{b} = sprintf('>%d', lo);
        else
            bin_labels{b} = sprintf('%d-%d', lo, hi);
        end
        fprintf('[%s: n=%d, JS=%.4f] ', bin_labels{b}, n_by_dist(b), js_by_dist(b));
    end
    fprintf('\n');
else
    bin_labels{1} = 'all';
    n_by_dist(1)   = sum(isfinite(js_vox));
    js_by_dist(1)  = nanmedian(js_vox);
    cos_by_dist(1) = nanmedian(cos_vox);
    rho_by_dist(1) = nanmedian(rho_vox);
end

% ---- 4. Summary ----
valid = isfinite(js_vox);
summary.median_js      = nanmedian(js_vox(valid));
summary.median_cosine  = nanmedian(cos_vox(valid));
summary.median_spearman= nanmedian(rho_vox(valid));
summary.n_valid        = sum(valid);

fprintf('  Overall: median JS=%.4f, cosine=%.4f, spearman=%.4f (n=%d)\n', ...
    summary.median_js, summary.median_cosine, summary.median_spearman, summary.n_valid);

% ---- output ----
val.js_voxelwise      = js_vox;
val.cosine_voxelwise  = cos_vox;
val.spearman_voxelwise= rho_vox;
val.dist_to_lesion    = dist_to_lesion;
val.bin_labels        = bin_labels;
val.js_by_dist        = js_by_dist;
val.cosine_by_dist    = cos_by_dist;
val.spearman_by_dist  = rho_by_dist;
val.n_by_dist         = n_by_dist;
val.summary           = summary;

end
