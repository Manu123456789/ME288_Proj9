function results = run_beta_validation(cfg, modeLookup, cases)
%RUN_BETA_VALIDATION Propagate and compare each independent truth case.

validateattributes(cfg.dt, {'numeric'}, {'real', 'finite', 'scalar', 'positive'}, ...
    mfilename, 'cfg.dt');
if ~isa(cfg.propagateCase, 'function_handle')
    error('BetaScaffold:Propagator', 'cfg.propagateCase must be a function handle.');
end
if ~iscell(modeLookup) || isempty(modeLookup) || ~isvector(modeLookup) || ...
        ~all(cellfun(@(lookup) isa(lookup, 'function_handle'), modeLookup))
    error('BetaScaffold:Tables', 'Supply one lookup function handle per mode.');
end
fields = {'id', 'R_0', 'V_0', 't_0', 't_final', 'truth'};
if ~isstruct(cases) || isempty(cases) || ~all(isfield(cases, fields))
    error('BetaScaffold:Cases', 'Supply cases with id, R_0, V_0, t_0, t_final, and truth.');
end
ids = {cases.id};
if ~all(cellfun(@(id) ischar(id) && isrow(id) && ~isempty(id), ids)) || ...
        numel(unique(ids)) ~= numel(ids)
    error('BetaScaffold:CaseIds', 'Case IDs must be unique nonempty character vectors.');
end
names = {cfg.channels.name};
if isempty(names) || ~all(cellfun(@isvarname, names)) || ...
        numel(unique(names)) ~= numel(names)
    error('BetaScaffold:Channels', 'Use unique MATLAB field names for validation channels.');
end

blank = struct('id', '', 't_0', [], 'truth', [], 'simulation', [], ...
    'switchTimes', [], 'comparison', [], 'status', '', 'message', '');
results.cases = repmat(blank, numel(cases), 1);
for k = 1:numel(cases)
    c = cases(k);
    item = blank;
    item.id = c.id;
    item.t_0 = c.t_0;
    item.truth = c.truth;
    try
        validateattributes(c.R_0, {'numeric'}, {'real', 'finite', 'vector', 'numel', 3});
        validateattributes(c.V_0, {'numeric'}, {'real', 'finite', 'vector', 'numel', 3});
        validateattributes(c.t_0, {'numeric'}, {'real', 'finite', 'scalar'});
        validateattributes(c.t_final, {'numeric'}, {'real', 'finite', 'scalar'});

        switchTimes = cfg.switchTimes;
        if isfield(c, 'switchTimes') && ~isempty(c.switchTimes)
            switchTimes = c.switchTimes;
        end
        validateattributes(switchTimes, {'numeric'}, {'real', 'finite', '2d'});
        if (~isempty(switchTimes) && ~isvector(switchTimes)) || ...
                numel(switchTimes) ~= numel(modeLookup) - 1 || ...
                any(diff(switchTimes(:)) <= 0)
            error('BetaScaffold:Schedule', ...
                'Case %s needs one increasing switch time per mode transition.', c.id);
        end
        item.switchTimes = switchTimes(:).';
        drag.switchTimes = item.switchTimes;
        drag.betaAt = @(t, query, segmentMode) ...
            estimate_beta(t, query, modeLookup, switchTimes, segmentMode);

        item.simulation = cfg.propagateCase(c.R_0(:), c.V_0(:), c.t_0, ...
            c.t_final, cfg.dt, drag);
        item.comparison = compare_beta_case(item.simulation, c, cfg.channels, ...
            item.switchTimes);
        item.status = item.comparison.status;
    catch problem
        if strcmp(problem.identifier, 'BetaScaffold:TODO')
            rethrow(problem);
        end
        item.status = 'ERROR';
        item.message = problem.message;
    end
    results.cases(k) = item;
end

compared = ~cellfun(@isempty, {results.cases.comparison});
results.metrics = struct([]);
if any(compared)
    comparisons = [results.cases(compared).comparison];
    results.metrics = vertcat(comparisons.metrics);
end
results.spread = summarize_cases(results.metrics, cfg.channels, numel(cases));
end

function spread = summarize_cases(metrics, channels, nCases)
% Each case contributes one value per metric and signal.
spread = struct([]);
measureNames = {'rms', 'peak', 'endpoint'};
for channelIndex = 1:numel(channels)
    channel = channels(channelIndex);
    selected = struct([]);
    if ~isempty(metrics)
        selected = metrics(strcmp({metrics.channel}, channel.name) & [metrics.complete]);
    end
    for k = 1:numel(measureNames)
        measure = measureNames{k};
        values = [];
        if ~isempty(selected)
            values = [selected.(measure)];
            values = values(isfinite(values));
        end
        row = struct('channel', channel.name, 'measure', measure, 'unit', channel.unit, ...
            'values', values(:), 'nValid', numel(values), 'nCases', nCases, ...
            'mean', NaN, 'median', NaN, 'minimum', NaN, 'maximum', NaN);
        if ~isempty(values)
            row.mean = mean(values);
            row.median = median(values);
            row.minimum = min(values);
            row.maximum = max(values);
        end
        spread(end+1, 1) = row; %#ok<AGROW>
    end
end
end
