function plot_distance_comparison(res, outpath)
% Bar plot: JS divergence by distance bin for all three methods

bins = res.val_corr.bin_labels;
Nb = numel(bins);

js = [res.val_pcdm.js_by_dist, ...
      res.val_uncorr.js_by_dist, ...
      res.val_corr.js_by_dist];

fig = figure('Position', [100 100 600 400], 'Visible', 'off');
b = bar(js);
b(1).FaceColor = [0.7 0.7 0.7];
b(2).FaceColor = [0.9 0.4 0.3];
b(3).FaceColor = [0.3 0.6 0.9];
set(gca, 'XTickLabel', bins, 'FontSize', 11);
xlabel('Distance from tumor (voxels)');
ylabel('Median JS divergence (lower = better)');
legend({'pCDM', 'iCDM-uncorrected', 'iCDM-corrected'}, 'Location', 'best');
title(sprintf('Structure-Function Agreement: %s', res.id));
grid on;

saveas(fig, fullfile(outpath, 'validation_distance_plot.png'));
close(fig);
fprintf('  Saved: %s\n', fullfile(outpath, 'validation_distance_plot.png'));
end


function aggregate_and_report(results, outdir)
% Aggregate JS divergence across patients per distance bin

Npat = numel(results);
if Npat == 0, return; end

Nb = numel(results{1}.val_corr.bin_labels);
bins = results{1}.val_corr.bin_labels;

js_pcdm   = NaN(Npat, Nb);
js_uncorr = NaN(Npat, Nb);
js_corr   = NaN(Npat, Nb);

for p = 1:Npat
    r = results{p};
    js_pcdm(p,:)   = r.val_pcdm.js_by_dist';
    js_uncorr(p,:) = r.val_uncorr.js_by_dist';
    js_corr(p,:)   = r.val_corr.js_by_dist';
end

fprintf('\n%-12s  %12s  %12s  %12s  %12s\n', ...
    'Dist bin', 'pCDM', 'Uncorrected', 'Corrected', 'Improvement');
fprintf('%s\n', repmat('-', 1, 65));

for b = 1:Nb
    m_p = nanmean(js_pcdm(:,b));
    m_u = nanmean(js_uncorr(:,b));
    m_c = nanmean(js_corr(:,b));
    improv = (m_u - m_c) / m_u * 100;
    fprintf('%-12s  %12.4f  %12.4f  %12.4f  %+11.1f%%\n', ...
        bins{b}, m_p, m_u, m_c, improv);
end

% Overall
fprintf('%s\n', repmat('-', 1, 65));
fprintf('%-12s  %12.4f  %12.4f  %12.4f  %+11.1f%%\n', ...
    'Overall', ...
    nanmean(js_pcdm(:)), ...
    nanmean(js_uncorr(:)), ...
    nanmean(js_corr(:)), ...
    (nanmean(js_uncorr(:)) - nanmean(js_corr(:))) / nanmean(js_uncorr(:)) * 100);

% Save summary table
T = table(bins(:), nanmean(js_pcdm,1)', nanmean(js_uncorr,1)', nanmean(js_corr,1)', ...
    'VariableNames', {'DistBin', 'pCDM', 'iCDM_uncorr', 'iCDM_corr'});
writetable(T, fullfile(outdir, 'aggregate_validation.csv'));
save(fullfile(outdir, 'aggregate_results.mat'), 'results', 'T');
fprintf('\nSaved: %s\n', fullfile(outdir, 'aggregate_validation.csv'));
end


function mkdir_if_needed(d)
if ~exist(d, 'dir'), mkdir(d); end
end
