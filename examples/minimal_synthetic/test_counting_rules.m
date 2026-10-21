function T = test_counting_rules(workdir)
% TEST_COUNTING_RULES  Executable specification of the streamline-to-count rules (docs/CONTRACT.md, Sec. 2).
%
% A 12 x 12 x 12 atlas (1 mm, voxel (i,j,k) at world (i-1, j-1, k-1) mm) holds three targets:
%   channel 1 : slab i = 1..2      (cortical target, "left")
%   channel 2 : slab i = 11..12    (cortical target, "right")
%   channel 3 : voxels i = 6, j = 6..7, k = 6..7   (subcortical structure inside the white matter)
% Seven hand-built streamlines exercise each rule; every resulting count is asserted.
if nargin < 1, workdir = fullfile(tempdir, 'icdm_counting_test'); end
if ~exist(workdir,'dir'), mkdir(workdir); end
dim = [12 12 12]; M = [eye(3) -ones(3,1); 0 0 0 1];
A = zeros(dim);
A(1:2,:,:) = 1; A(11:12,:,:) = 2; A(6,6:7,6:7) = 3;
V = struct('fname',fullfile(workdir,'atlas.nii'),'dim',dim,'dt',[spm_type('uint8') 0],'mat',M,'pinfo',[1;0;0],'n',[1 1]);
spm_write_vol(V, A);

w = @(i,j,k) [i-1; j-1; k-1];                     % voxel -> world (mm)
line = @(i0,i1,j,k) cell2mat(arrayfun(@(i) w(i,j,k), i0:0.5:i1, 'uni', 0));
sl = {};
sl{1} = line(2,11,3,3);                           % A: endpoint to endpoint            -> {1,2}
sl{2} = line(4,10,4,4);                           % B: both ends one voxel short       -> snapped, {1,2}
sl{3} = line(2,11,6,6);                           % C: passes subcortical channel 3    -> {1,2,3}
sl{4} = line(5,8,9,9);                            % D: never reaches a label           -> dropped
sl{5} = [line(2,5,10,3) w(5,11,3) w(5,10,3) line(6,11,10,3)];   % E: re-enters voxel (5,10,3) -> {1,2}
sl{6} = line(2,7,8,8);                            % F: one end unassigned              -> {1}
sl{7} = [w(-5,2,10) w(-2,2,10) line(1,11,2,10)];  % G: starts outside the grid         -> {1,2}
write_tck(fullfile(workdir,'rules.tck'), sl);

C = double(icdm_tracks_to_counts(fullfile(workdir,'rules.tck'), V.fname, fullfile(workdir,'counts.nii')));
n = @(i,j,k) squeeze(C(i,j,k,:))';

% A: every visited voxel counts once for channels 1 and 2
assert(isequal(n(2,3,3),[1 1 0]) && isequal(n(7,3,3),[1 1 0]) && isequal(n(11,3,3),[1 1 0]));
% B: unlabelled ends are extended along the terminal direction; the labelled voxels reached are added,
%    the unlabelled gap voxel is not
assert(isequal(n(2,4,4),[1 1 0]) && isequal(n(3,4,4),[0 0 0]) && isequal(n(4,4,4),[1 1 0]));
assert(isequal(n(11,4,4),[1 1 0]) && isequal(n(12,4,4),[1 1 0]));
% C: a labelled structure on the path is a target; it is counted in every visited voxel, itself included
assert(isequal(n(4,6,6),[1 1 1]) && isequal(n(6,6,6),[1 1 1]));
% D: a streamline with no label contributes nothing
assert(all(C(:,9,9,:) == 0, 'all'));
% E: consecutive repeats collapse; re-entering a voxel after leaving it counts again
assert(isequal(n(4,10,3),[1 1 0]) && isequal(n(5,10,3),[2 2 0]) && isequal(n(5,11,3),[1 1 0]));
% F: one assigned end -> counted for that target only
assert(isequal(n(3,8,8),[1 0 0]) && isequal(n(7,8,8),[1 0 0]));
% G: vertices outside the grid are dropped
assert(isequal(n(1,2,10),[1 1 0]) && isequal(n(11,2,10),[1 1 0]));
% each target is counted once per streamline: channel totals over the grid
T.total_per_channel = squeeze(sum(sum(sum(C,1),2),3))';
T.N_max = max(sum(C,4), [], 'all');
fprintf('[test_counting_rules] all assertions passed; totals per channel = %s\n', mat2str(T.total_per_channel));
end
