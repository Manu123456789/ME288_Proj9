function out = plot_drag_cd_map(time, R, V, dragFcn, reDivisor, markerStride)
%PLOT_DRAG_CD_MAP Plot wind-axis CD and the queries along a trajectory.
%   out = plot_drag_cd_map(time, R, V)
%   out = plot_drag_cd_map(time, R, V, @Drag, 1e6)
%   out = plot_drag_cd_map(time, R, V, @Drag, 1e6, 5)
%
% R and V are N-by-3 ECI histories. Pass rows in propagation order.
% Use the Drag function with the optional third output described in README.
% reDivisor changes the display axis only. The default is 1 (table scale).
% markerStride shows every Nth valid query as a red dot. The default is 1.
% The line, arrows, and returned arrays use all saved states.
% out.Mach, out.Re, and out.CD contain the evaluated values for each row.
% A state with zero relative airspeed has NaN entries and breaks the line.

if nargin < 4 || isempty(dragFcn)
    dragFcn = @Drag;
end
if nargin < 5 || isempty(reDivisor)
    reDivisor = 1;
end
if nargin < 6 || isempty(markerStride)
    markerStride = 1;
end
validateattributes(dragFcn, {'function_handle'}, {'scalar'});
validateattributes(reDivisor, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(markerStride, {'numeric'}, {'scalar','real','finite','integer','positive'});
reDivisor = double(reDivisor);
markerStride = double(markerStride);
validateattributes(time, {'numeric'}, {'vector','real','finite','nonempty'});
time = double(time(:));
n = numel(time);
validateattributes(R, {'numeric'}, {'size',[n 3],'real','finite'});
validateattributes(V, {'numeric'}, {'size',[n 3],'real','finite'});

% Replay the saved states through the same Drag used by the propagator.
mach = nan(n,1);
reynolds = nan(n,1);
cd = nan(n,1);
for k = 1:n
    if k == 1
        [~, aero, map] = dragFcn(R(k,:), V(k,:));
    else
        [~, aero] = dragFcn(R(k,:), V(k,:));
    end
    if ~isempty(aero)
        mach(k) = aero.Mach;
        reynolds(k) = aero.Re;
        cd(k) = aero.CD;
    end
end

rePlot = reynolds / reDivisor;
valid = isfinite(mach) & isfinite(reynolds) & isfinite(cd);
ink = [0.11, 0.19, 0.26];
red = [0.78, 0.11, 0.17];
fig = figure('Name', 'Drag coefficient query path', 'Color', 'w', ...
    'Position', [110 100 970 690]);
ax = axes('Parent', fig, 'Color', [0.96, 0.97, 0.98]);
hold(ax, 'on');
drawCoefficientSurface(ax, map, reDivisor);

if any(valid)
    % Keep NaNs in the line so it breaks across states with no lookup.
    plot3(ax, mach, rePlot, cd, '-', 'Color', 'w', 'LineWidth', 3.2, ...
        'HandleVisibility', 'off');
    path = plot3(ax, mach, rePlot, cd, '-', ...
        'Color', [0.49, 0.09, 0.13], 'LineWidth', 1.65, ...
        'DisplayName', 'Propagation order');
    dotRows = find(valid);
    dotRows = unique([dotRows(1:markerStride:end); dotRows(end)]);
    dots = scatter3(ax, mach(dotRows), rePlot(dotRows), cd(dotRows), 26, 'filled', ...
        'MarkerFaceColor', red, 'MarkerEdgeColor', 'w', 'LineWidth', 0.7, ...
        'DisplayName', 'Saved-state drag queries');
    dots.DataTipTemplate.DataTipRows = [ ...
        datatipTextRow('Time (s)', time(dotRows)), ...
        datatipTextRow('Mach', mach(dotRows)), ...
        datatipTextRow('Reynolds (table scale)', reynolds(dotRows)), ...
        datatipTextRow('C_D', cd(dotRows))];
    first = dotRows(1);
    last = dotRows(end);
    if first ~= last && (mach(first) ~= mach(last) || rePlot(first) ~= rePlot(last))
        text(ax, mach(first), rePlot(first), cd(first), 'First', ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
            'Color', ink, 'BackgroundColor', 'w', 'Margin', 2, 'FontWeight', 'bold', ...
            'Clipping', 'on', 'HitTest', 'off', 'PickableParts', 'none');
        text(ax, mach(last), rePlot(last), cd(last), 'Last', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
            'Color', ink, 'BackgroundColor', 'w', 'Margin', 2, 'FontWeight', 'bold', ...
            'Clipping', 'on', 'HitTest', 'off', 'PickableParts', 'none');
    end
    legend(ax, [dots path], 'Location', 'northeast');
end

view(ax, 2);
set(ax, 'SortMethod', 'childorder', 'Layer', 'top', 'FontSize', 11, ...
    'XColor', ink, 'YColor', ink, 'GridAlpha', 0.16);
grid(ax, 'on');
box(ax, 'on');
colormap(ax, parula(256));
key = colorbar(ax);
key.Label.String = 'Wind-axis drag coefficient C_D';
xlabel(ax, 'Mach number');
if reDivisor == 1
    ylabel(ax, 'Reynolds number (table scale)');
elseif reDivisor == 1e6
    ylabel(ax, 'Table Reynolds number (millions)');
else
    ylabel(ax, sprintf('Table Reynolds number / %g', reDivisor));
end
zlabel(ax, 'C_D');
title(ax, sprintf('Drag coefficient at %.3g deg angle of attack', map.alphaDeg));
drawDirectionArrows(ax, mach, rePlot, cd, valid, red);

out = struct('figure', fig, 'axes', ax, 'time', time, ...
    'Mach', mach, 'Re', reynolds, 'CD', cd);
end

function drawCoefficientSurface(ax, map, reDivisor)
% Partition the surface at every Reynolds breakpoint of an enabled table.
tables = map.tables(~cellfun(@isempty, map.tables));
grids = cellfun(@(table) table.Re(:), tables, 'UniformOutput', false);
first = max(cellfun(@(grid) grid(1), grids));
last = min(cellfun(@(grid) grid(end), grids));
if first > last
    error('DragCDMap:NoSharedReynolds', ...
        'The enabled tables have no common Reynolds range.');
end
edges = unique(vertcat(grids{:}));
edges = edges(edges >= first & edges <= last);

for k = 1:numel(edges)-1
    range = edges(k:k+1);
    mach = sharedMachBreakpoints(tables, range);
    if isempty(mach)
        continue
    end
    mach = unique([mach; mach(1:end-1) + diff(mach)/2]);
    re = linspace(range(1), range(2), 4);
    [M, Re] = meshgrid(mach, re);
    CD = sampleCD(map, M, Re);
    if numel(mach) == 1
        keep = isfinite(CD);
        scatter3(ax, M(keep), Re(keep)/reDivisor, CD(keep), 14, CD(keep), 'filled', ...
            'HandleVisibility', 'off');
    else
        surface(ax, M, Re/reDivisor, CD, CD, 'FaceColor', 'interp', ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

% Exact Reynolds breakpoints can cover Mach values beyond adjacent strips.
for k = 1:numel(edges)
    re = edges(k);
    mach = sharedMachBreakpoints(tables, re);
    if isempty(mach)
        continue
    end
    re = repmat(re, size(mach));
    cd = sampleCD(map, mach, re);
    keep = isfinite(cd);
    scatter3(ax, mach(keep), re(keep)/reDivisor, cd(keep), 9, cd(keep), 'filled', ...
        'MarkerEdgeColor', [0.18, 0.24, 0.30], 'MarkerEdgeAlpha', 0.25, ...
        'HandleVisibility', 'off');
end
end

function breakpoints = sharedMachBreakpoints(tables, reRange)
% Find the Mach coverage shared by the slices used over this Re interval.
lower = -Inf;
upper = Inf;
breakpoints = [];
for k = 1:numel(tables)
    table = tables{k};
    first = find(table.Re <= reRange(1), 1, 'last');
    last = find(table.Re >= reRange(end), 1, 'first');
    for j = first:last
        curve = table.curves{j};
        mach = curve.GridVectors{1}(:);
        lower = max(lower, mach(1));
        upper = min(upper, mach(end));
        breakpoints = [breakpoints; mach]; %#ok<AGROW>
    end
end
breakpoints = unique(breakpoints(breakpoints >= lower & breakpoints <= upper));
end

function CD = sampleCD(map, Mach, Re)
CD = arrayfun(map.evaluate, Mach, Re);
CD(~isfinite(CD) | CD < 0) = NaN;
end

function drawDirectionArrows(ax, mach, re, cd, valid, color)
% Use distance in normalized axis coordinates to space the arrowheads.
xLimits = xlim(ax);
yLimits = ylim(ax);
span = [diff(xLimits), diff(yLimits)];
origin = [xLimits(1), yLimits(1)];
xy = ([mach, re] - origin) ./ span;
step = diff(xy, 1, 1);
lengths = hypot(step(:,1), step(:,2));
segments = find(valid(1:end-1) & valid(2:end) & lengths > 0);
if isempty(segments)
    return
end

% Keep the arrowheads from expanding the data limits.
xlim(ax, xLimits);
ylim(ax, yLimits);
distance = cumsum(lengths(segments));
total = distance(end);
count = min(5, max(1, floor(total / 0.10)));
spacing = total / (count + 1);
headLength = min(0.023, 0.32 * spacing);
for target = spacing * (1:count)
    j = find(distance >= target, 1, 'first');
    k = segments(j);
    direction = step(k,:) / lengths(k);
    across = [-direction(2), direction(1)];
    fraction = (target - (distance(j) - lengths(k))) / lengths(k);
    tip = xy(k,:) + fraction * step(k,:);
    base = tip - headLength * direction;
    points = [tip; base + 0.43 * headLength * across; ...
                   base - 0.43 * headLength * across];
    points = points .* span + origin;
    height = cd(k) + fraction * (cd(k+1) - cd(k));
    patch(ax, 'XData', points(:,1), 'YData', points(:,2), ...
        'ZData', repmat(height, 3, 1), 'FaceColor', color, ...
        'EdgeColor', 'w', 'LineWidth', 0.65, 'HandleVisibility', 'off', ...
        'HitTest', 'off', 'PickableParts', 'none');
end
end
