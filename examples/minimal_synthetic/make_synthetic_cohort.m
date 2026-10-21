function info = make_synthetic_cohort(outdir, S, nstream, seed, age)
% MAKE_SYNTHETIC_COHORT  Small synthetic tractograms + atlas for the minimal end-to-end example.
%
%   info = make_synthetic_cohort(outdir, S, nstream, seed, age)
%   age (optional, S x 1, years on an 8-12 scale) sets the covariate that drives the synthetic effect
%
% Writes, in outdir:
%   atlas.nii      30 x 30 x 20 grid, 2 mm, RAS; labels 1-6 = cortical targets on the six faces,
%                  label 7 = a subcortical structure in the white matter (counted, but not a model target)
%   template.nii   the analysis (template) grid; here identical to the atlas grid
%   sub-XX/tracts.tck              MRtrix tractogram in world mm (nstream streamlines)
%   sub-XX/dartel/y_fs_t1.nii      identity SPM deformation (subject space == template space)
%   participants.tsv               participant_id, age
% Streamline target pairs are drawn from a subject-specific mixture whose weight on the left-right pair
% increases with age, so the cohort carries a known compositional age effect.
if nargin < 2, S = 8; end
if nargin < 3, nstream = 3000; end
if nargin < 4, seed = 1; end
rng(seed);
if ~exist(outdir,'dir'), mkdir(outdir); end
dim = [30 30 20]; vox = 2;
M = [vox 0 0 -vox*(dim(1)+1)/2; 0 vox 0 -vox*(dim(2)+1)/2; 0 0 vox -vox*(dim(3)+1)/2; 0 0 0 1];

A = zeros(dim);
A(:,:,1:2) = 6; A(:,:,end-1:end) = 5;
A(:,1:3,:) = 3; A(:,end-2:end,:) = 4;
A(1:3,:,:) = 1; A(end-2:end,:,:) = 2;
A(14:17,14:17,9:12) = 7;
Vt = struct('fname',fullfile(outdir,'atlas.nii'),'dim',dim,'dt',[spm_type('uint8') 0],'mat',M,'pinfo',[1;0;0],'n',[1 1]);
spm_write_vol(Vt, A);
Vt.fname = fullfile(outdir,'template.nii'); Vt.dt = [spm_type('float32') 0];
spm_write_vol(Vt, ones(dim));

pairs = nchoosek(1:6, 2);                       % 15 cortical target pairs
base = ones(1, size(pairs,1)); base(1) = 4;     % pair (1,2) is the most frequent
if nargin < 5 || isempty(age), age = round(linspace(8, 12, S)' + 0.3*randn(S,1), 1); end
ctr = (dim + 1) / 2;
fid = fopen(fullfile(outdir,'participants.tsv'),'w'); fprintf(fid,'participant_id\tage\n');
for s = 1:S
    id = sprintf('sub-%02d', s);
    sd = fullfile(outdir, id); if ~exist(fullfile(sd,'dartel'),'dir'), mkdir(fullfile(sd,'dartel')); end
    w = base .* exp(0.3*randn(size(base))); w(1) = w(1) * (1 + 0.15*(age(s)-10));
    w = w / sum(w);
    sl = cell(nstream, 1);
    for t = 1:nstream
        u = rand;
        if u < 0.03                                 % short fragment in white matter: no cortical target
            p0 = ctr + [randi([-5 5]) randi([-5 5]) randi([-3 3])];
            sl{t} = curve(p0, p0 + [3 0 0], p0 + [6 0 0]);
            continue
        end
        pr = pairs(find(rand < cumsum(w), 1), :);
        p0 = face_point(pr(1), dim); p1 = face_point(pr(2), dim);
        if u < 0.13                                 % routed through the subcortical structure
            c = [15.5 15.5 10.5];
        else
            c = ctr + [2 2 1].*randn(1,3);
        end
        sl{t} = curve(p0, c, p1);
    end
    sl = cellfun(@(v) M(1:3,:) * [v; ones(1,size(v,2))], sl, 'uni', 0);   % voxel -> world mm
    write_tck(fullfile(sd,'tracts.tck'), sl);
    make_identity_deformation(fullfile(outdir,'template.nii'), fullfile(sd,'dartel','y_fs_t1.nii'));
    fprintf(fid, '%s\t%.1f\n', id, age(s));
end
fclose(fid);
info = struct('outdir',outdir,'S',S,'age',age,'dim',dim,'mat',M,'targets',1:6);
end

function p = face_point(lab, dim)
% a point inside the target slab of label lab
r = @(n) 4 + rand*(n - 7);
switch lab
    case 1, p = [2 r(dim(2)) r(dim(3))];
    case 2, p = [dim(1)-1 r(dim(2)) r(dim(3))];
    case 3, p = [r(dim(1)) 2 r(dim(3))];
    case 4, p = [r(dim(1)) dim(2)-1 r(dim(3))];
    case 5, p = [r(dim(1)) r(dim(2)) dim(3)-0.5];
    case 6, p = [r(dim(1)) r(dim(2)) 1.5];
end
end

function v = curve(p0, c, p1)
% quadratic Bezier from p0 through control point c to p1, sampled every ~0.5 voxel
L = norm(c - p0) + norm(p1 - c); n = max(2, ceil(2*L));
t = linspace(0, 1, n);
v = ((1-t).^2)'*p0 + (2*(1-t).*t)'*c + (t.^2)'*p1;
v = v';
end
