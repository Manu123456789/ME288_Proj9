function comparison = compare_beta_case(sim, c, channels, switchTimes)
%COMPARE_BETA_CASE Compare simulation output at the original truth timestamps.

[tSim, simOrder] = ordered_time(sim);
[tRef, refOrder] = ordered_time(c.truth);
window = sort([double(c.t_0), double(c.t_final)]);
inWindow = tRef >= window(1) & tRef <= window(2);
queryTime = tRef(inWindow);
covered = tSim(1) <= window(1) && tSim(end) >= window(2) && ...
    tRef(1) <= window(1) && tRef(end) >= window(2);

comparison.time = queryTime;
comparison.errors = struct();
blank = struct('id', c.id, 'channel', '', 'unit', '', 'required', false, ...
    'nSamples', numel(queryTime), 'nMatched', 0, 'complete', false, ...
    'rms', NaN, 'peak', NaN, 'endpoint', NaN, 'peakLimit', NaN, ...
    'margin', NaN, 'status', 'MISSING');
metrics = repmat(blank, numel(channels), 1);
for k = 1:numel(channels)
    channel = channels(k);
    name = channel.name;
    row = blank;
    row.channel = name;
    row.unit = channel.unit;
    row.required = channel.required || ~isempty(channel.peakLimit);
    if ~isempty(channel.peakLimit)
        validateattributes(channel.peakLimit, {'numeric'}, ...
            {'real', 'finite', 'scalar', 'nonnegative'});
        row.peakLimit = double(channel.peakLimit);
    end

    if isfield(sim.signals, name) && isfield(c.truth.signals, name)
        if isempty(channel.unit)
            row.status = 'UNITS_UNSET';
        else
            predicted = read_signal(sim, channel, simOrder);
            reference = read_signal(c.truth, channel, refOrder(inWindow));
            breaks = [];
            if channel.splitAtSwitch
                breaks = switchTimes;
            end
            delta = sample_history(tSim, predicted, queryTime, breaks) - reference;
            valid = all(isfinite(delta), 2);
            delta(~valid, :) = NaN;
            magnitude = vecnorm(delta, 2, 2);
            comparison.errors.(name) = struct('components', delta, 'magnitude', magnitude);
            row.nMatched = nnz(valid);
            row.complete = covered && ~isempty(queryTime) && all(valid);
            if any(valid)
                row.rms = norm(magnitude(valid)) / sqrt(row.nMatched);
                row.peak = max(magnitude(valid));
            end
            endpoint = queryTime == c.t_final & valid;
            if any(endpoint)
                row.endpoint = magnitude(endpoint);
            end
            row.margin = row.peakLimit - row.peak;
            if row.margin < 0
                row.status = 'FAIL';
            elseif ~row.complete
                row.status = 'INCOMPLETE';
            elseif isnan(row.peakLimit)
                row.status = 'UNASSESSED';
            else
                row.status = 'PASS';
            end
        end
    end
    metrics(k) = row;
end
comparison.metrics = metrics;
required = metrics([metrics.required]);
if any(strcmp({required.status}, 'FAIL'))
    comparison.status = 'FAIL';
elseif any(~[required.complete])
    comparison.status = 'INCOMPLETE';
elseif isempty(required) || any(strcmp({required.status}, 'UNASSESSED'))
    comparison.status = 'UNASSESSED';
else
    comparison.status = 'PASS';
end
end

function [t, order] = ordered_time(history)
validateattributes(history.time, {'numeric'}, {'real', 'finite', 'vector', 'nonempty'});
if ~isstruct(history.signals) || ~isscalar(history.signals)
    error('BetaScaffold:Signals', 'history.signals must be a scalar struct.');
end
[t, order] = sort(double(history.time(:)));
if any(diff(t) == 0)
    error('BetaScaffold:Time', 'History timestamps must be unique.');
end
end

function values = read_signal(history, channel, order)
values = history.signals.(channel.name);
validateattributes(values, {'numeric'}, ...
    {'real', 'size', [numel(history.time), channel.width]}, mfilename, channel.name);
values = double(values(order, :));
values(~isfinite(values)) = NaN;
end

function aligned = sample_history(t, values, query, breaks)
aligned = nan(numel(query), size(values, 2));
edges = [-Inf, breaks(:).', Inf];
for k = 1:numel(edges)-1
    source = t >= edges(k) & t < edges(k+1);
    target = query >= edges(k) & query < edges(k+1);
    if nnz(source) > 1 && any(target)
        aligned(target, :) = interp1(t(source), values(source, :), ...
            query(target), 'linear', NaN);
    end
end
% Exact samples also cover segments with a single available simulation row.
[exact, index] = ismember(query, t);
aligned(exact, :) = values(index(exact), :);
end
