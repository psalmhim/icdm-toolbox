function check_outputs(summary_file, expected_file, rtol)
% CHECK_OUTPUTS  Compare a run's summary.json with the reference values in expected_outputs.json.
% Integer fields (voxel counts, N_total) must match exactly; medians within relative tolerance rtol (1e-4).
if nargin < 3, rtol = 1e-4; end
a = jsondecode(fileread(summary_file)); b = jsondecode(fileread(expected_file));
assert(a.n_group_voxels == b.n_group_voxels, 'group mask size differs: %d vs %d', a.n_group_voxels, b.n_group_voxels);
f_exact = {'n_native_valid','n_template_valid','N_total'};
f_tol = {'kappa_data_median','kappa_post_median','r_data_median'};
for s = 1:numel(b.subjects)
    x = a.subjects(s); y = b.subjects(s);
    assert(strcmp(x.id, y.id));
    for q = f_exact, assert(x.(q{1}) == y.(q{1}), '%s %s: %g vs %g', x.id, q{1}, x.(q{1}), y.(q{1})); end
    for q = f_tol
        assert(abs(x.(q{1}) - y.(q{1})) <= rtol*max(1,abs(y.(q{1}))), '%s %s: %g vs %g', x.id, q{1}, x.(q{1}), y.(q{1}));
    end
end
fprintf('[check_outputs] %d subjects match the reference values\n', numel(b.subjects));
end
