function gopts = icdm_evaluate_group_mask(subjects, min_coverage_frac, outfile)
% =========================================================================
% icdm_evaluate_group_mask.m
%
% 입력:
%   subjects          : struct array (icdm_compose_subject 결과)
%   min_coverage_frac : 예) 0.3
%   outfile           : output file for saving results
% 필요:
%   각 subj.mni_features_mat 에서
%       dim_mni, idx_mni_mask를 읽어 coverage 계산
%
% 출력:
%   gopts.dim_mni      : [1x3] MNI dimension
%   gopts.coverage_mni : logical vector (prod(dim_mni) x 1)
%   gopts.idx_mni      : find(coverage_mni)
% =========================================================================
if nargin<3, outfile = ''; end
if nargin<2, min_coverage_frac = 0.3; end
S = numel(subjects);
if S == 0
    error('No subjects given.');
end

if exist(outfile,'file')
    fprintf('[GroupMask] Loading existing group mask file: %s\n', outfile);
    L = load(outfile, 'dim_mni','idx_mni');
    gopts.dim_mni = L.dim_mni;
    gopts.idx_mni = L.idx_mni;
    return;
end

fprintf('[GroupMask] Evaluating group MNI coverage over %d subjects...\n', S);
% 각 subject의 MNI mask index를 coverage에 더함
cov_count=[];
for s=1:S
    if rem(s,10)==1
        fprintf('  Processing subject %d/%d: %s\n', s, S, subjects(s).id);
    end
    info = load(subjects(s).datafile);
    M  = zeros(info.dim_native,"single");
    M(info.idx_native) = 1;
    % native -> MNI (nearest neighbor)
    Yout = icdm_warp_4d(M, info.V, subjects(s).dartel_flow, subjects(s).template, +1, 0);
    dim_mni = size(Yout);
    if numel(dim_mni) < 3
        dim_mni(3) = 1;
    end
    idx_mask = find(Yout > 0.5);
    if s==1
        Nvox_mni = prod(dim_mni);
        cov_count = zeros(Nvox_mni,1,'uint16');
    end
    cov_count(idx_mask) = cov_count(idx_mask) + 1;
end

idx_mni = uint32(find(cov_count / S >= min_coverage_frac));

gopts.dim_mni      = dim_mni;
gopts.idx_mni      = uint32(idx_mni);

fprintf('[GroupMask] dim=[%d %d %d], coverage voxels=%d (%.2f%%)\n', ...
    dim_mni(1),dim_mni(2),dim_mni(3), ...
    numel(idx_mni), 100*numel(idx_mni)/Nvox_mni);

if ~isempty(outfile)
    save(outfile,'idx_mni','dim_mni', '-v7.3');
    fprintf('[GroupMask] Saved group mask info to %s\n', outfile);
end
end

