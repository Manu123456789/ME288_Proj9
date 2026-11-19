function out = fit_loft_drag(runs, dt, scaleGrid, rTol, vTol, holdoutLofts)
% R (N-by-3 m) and V (N-by-3 m/s) are ECI. Time is in seconds.
% runs(1) = struct('label','A','loftDeg',15,'time',tA,'R',RA,'V',VA);
% Add one struct per trajectory. Repeated lofts share a fitted scale.
% In Drag.m use: inverseBeta = drag_study_scale() * model.K_norm * CD;
% holdoutLofts is optional. rTol (m) and vTol (m/s) weight the errors.
% Example: out = fit_loft_drag(runs,dt,0.90:0.025:1.10,100,1);
% out.fit gives one scale per training loft. out.cases compares run scores.
if nargin < 6, holdoutLofts = []; end
runs = runs(:);
loft = [runs.loftDeg].';
train = ~ismember(loft, holdoutLofts);
levels = unique(loft(train));
g = unique(scaleGrid(:)).';
assert(any(train) && all(isfinite(loft)), 'Supply at least one training loft and finite loft angles.');
assert(numel(g) >= 2 && all(isfinite(g) & g > 0), 'Supply at least two distinct positive finite scales.');
validateattributes(rTol, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(vTol, {'numeric'}, {'scalar','real','finite','positive'});
if g(1) < 1 && g(end) > 1, g = unique([g 1]); end

old = drag_study_scale();
restore = onCleanup(@() drag_study_scale(old)); %#ok<NASGU>
score = @(s,cases) caseLoss(s, cases, dt, rTol, vTol);
training = runs(train);
groups = [true(numel(training),1), loft(train) == levels.']; % pooled, then each loft
gridLoss = zeros(numel(training), numel(g));
for j = 1:numel(g)
    gridLoss(:,j) = score(g(j), training);
end

best = zeros(size(groups,2),1);
profiles = zeros(size(groups,2),numel(g));
for k = 1:numel(best)
    use = groups(:,k);
    y = mean(gridLoss(use,:),1);
    profiles(k,:) = y;
    tied = find(y == min(y));
    [~,nearest] = min(abs(g(tied) - 1));
    j = tied(nearest);
    best(k) = g(j);
    if all(y == y(1)), continue, end
    f = @(s) mean(score(s,training(use)));
    [s,loss] = fminbnd(f, g(max(1,j-1)), g(min(numel(g),j+1)));
    if loss < y(j), best(k) = s; end
end
out.globalScale = best(1);
out.fit = table(levels,best(2:end),'VariableNames',{'LoftDeg','Scale'});
out.grid = g;
out.sweep = profiles(2:end,:); % rows match out.fit, columns match out.grid

% Per-run comparison. A held-out loft outside the training range gets NaN.
sLoft = nan(size(loft));
if numel(levels) > 1
    sLoft = interp1(levels, best(2:end), loft, 'linear', NaN);
else
    sLoft(loft == levels(1)) = best(2);
end
Jloft = nan(size(loft));
for i = find(isfinite(sLoft)).'
    Jloft(i) = score(sLoft(i),runs(i));
end
out.cases = table(string({runs.label}).',loft,~train,sLoft,score(1,runs),score(best(1),runs),Jloft, ...
    'VariableNames',{'Label','LoftDeg','HeldOut','LoftScale','Jbase','Jglobal','Jloft'});
end

function J = caseLoss(scale, cases, dt, rTol, vTol)
drag_study_scale(scale);
J = zeros(numel(cases),1);
for i = 1:numel(cases)
    c = cases(i);
    t = c.time(:);
    assert(numel(t) >= 2 && (all(diff(t)>0) || all(diff(t)<0)) && ...
        isequal(size(c.R),size(c.V),[numel(t),3]), 'Check time order and N-by-3 R,V histories.');
    [R,V,tp] = propagate_traj(c.R(1,:),c.V(1,:),t(1),t(end),dt);
    tp = tp(:);
    tol = 32*eps(max([1;abs(t([1 end]))]));
    if abs(tp(1)-t(1)) <= tol, tp(1) = t(1); end
    if abs(tp(end)-t(end)) <= tol, tp(end) = t(end); end
    e = interp1(tp,[R V],t,'linear') - [c.R c.V];
    z = sum(e(:,1:3).^2,2)/rTol^2 + sum(e(:,4:6).^2,2)/vTol^2;
    J(i) = trapz(t,z)/(t(end)-t(1));
    assert(isfinite(J(i)), 'Check state values and propagated time coverage for %s.',char(c.label));
end
end
