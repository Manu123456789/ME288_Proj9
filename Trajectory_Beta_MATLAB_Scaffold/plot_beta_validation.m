function plot_beta_validation(results, channel)
%PLOT_BETA_VALIDATION ECI overlays, error histories, and case RMS distribution.
if nargin < 2
    channel = 'R';
end
rmsRows = results.spread(strcmp({results.spread.channel}, channel) & ...
    strcmp({results.spread.measure}, 'rms'));
if isempty(rmsRows)
    error('BetaScaffold:PlotChannel', 'Choose a configured validation channel.');
end

figure('Name', 'Trajectory beta validation');
layout = tiledlayout(1, 3);
pathAxes = nexttile(layout);
hold(pathAxes, 'on');
errorAxes = nexttile(layout);
hold(errorAxes, 'on');
distributionAxes = nexttile(layout);
pair = [];
for k = 1:numel(results.cases)
    c = results.cases(k);
    if isempty(c.comparison)
        continue;
    end
    if isfield(c.comparison.errors, 'R')
        r = c.truth.signals.R;
        p = c.simulation.signals.R;
        pair = [plot3(pathAxes, r(:,1), r(:,2), r(:,3), 'Color', [0.5 0.5 0.5]), ...
                plot3(pathAxes, p(:,1), p(:,2), p(:,3), 'b--')];
    end
    if isfield(c.comparison.errors, channel)
        elapsed = abs(c.comparison.time - double(c.t_0));
        plot(errorAxes, elapsed, c.comparison.errors.(channel).magnitude);
    end
end
title(pathAxes, 'ECI trajectory overlays');
xlabel(pathAxes, 'X (m)'); ylabel(pathAxes, 'Y (m)'); zlabel(pathAxes, 'Z (m)');
axis(pathAxes, 'equal'); view(pathAxes, 3); grid(pathAxes, 'on');
if ~isempty(pair)
    legend(pathAxes, pair, {'Truth', 'Model'}, 'Location', 'best');
end
title(errorAxes, 'Error histories by case');
xlabel(errorAxes, 'Elapsed time (s)');
ylabel(errorAxes, sprintf('%s error magnitude (%s)', channel, rmsRows.unit));
grid(errorAxes, 'on');

values = sort(rmsRows.values);
if ~isempty(values)
    stairs(distributionAxes, [values(1); values], (0:numel(values)).' / numel(values));
end
title(distributionAxes, sprintf('Case RMS spread (%d/%d)', rmsRows.nValid, rmsRows.nCases));
xlabel(distributionAxes, sprintf('%s RMS error (%s)', channel, rmsRows.unit));
ylabel(distributionAxes, 'Empirical cumulative fraction');
ylim(distributionAxes, [0 1]); grid(distributionAxes, 'on');
end
