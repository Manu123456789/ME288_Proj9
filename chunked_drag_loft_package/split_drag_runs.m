function chunks = split_drag_runs(runs,chunkSeconds)
% Split truth histories into independently initialized calibration windows.
% Keep the original absolute time and loft, including exact window edges.
if ~isscalar(chunkSeconds) || ~isfinite(chunkSeconds) || chunkSeconds<=0
    error('chunkSeconds must be positive and finite.');
end
if isempty(runs), error('At least one truth run is required.'); end
chunks = struct([]);
for j = 1:numel(runs)
    r = runs(j);
    t = r.time(:);
    if numel(t)<2 || any(~isfinite(t)) || ...
            ~(all(diff(t)>0) || all(diff(t)<0)) || ...
            ~isequal(size(r.R),[numel(t),3]) || ...
            ~isequal(size(r.V),[numel(t),3]) || ...
            any(~isfinite([r.R(:);r.V(:)]))
        error('Run %d needs monotonic time and finite Nx3 R and V.',j);
    end
    direction = sign(t(end)-t(1));
    duration = abs(t(end)-t(1));
    edges = t(1) + direction*(0:floor(duration/chunkSeconds))*chunkSeconds;
    tol = min(duration/4,16*eps(max(1,max(abs(t)))));
    if numel(edges)>1 && abs(edges(end)-t(end))<=tol
        edges(end) = t(end);
    else
        edges(end+1) = t(end);
    end
    [ts,ix] = sort(t);
    for k = 1:numel(edges)-1
        a = edges(k); b = edges(k+1);
        inside = direction*(t-a)>0 & direction*(b-t)>0;
        tc = [a;t(inside);b];
        % >>> NEW: both initial states come from truth, never previous output.
        c.label = sprintf('%s | chunk %d',char(r.label),k);
        c.loftDeg = r.loftDeg;
        c.time = tc;
        c.R = interp1(ts,r.R(ix,:),tc,'linear');
        c.V = interp1(ts,r.V(ix,:),tc,'linear');
        c.parentCase = j;
        c.chunkIndex = k;
        c.interval = [a b];
        chunks(end+1) = c; %#ok<AGROW>
    end
end
end
