function out = fit_loft_drag_time(runs, dt, nKnots)
%FIT_LOFT_DRAG_TIME Fit time-varying drag scales for each truth trajectory.
% runs(j): label, loftDeg, time (Nx1), R (Nx3), V (Nx3).
% Default: optimize at every truth time. Set nKnots=25 for a faster fit.
% Requires Optimization Toolbox (lsqnonlin).

opts = optimoptions('lsqnonlin','Display','off');
cleanup = onCleanup(@() drag_study_scale('reset'));
out = struct([]);

for j = 1:numel(runs)
    run = runs(j);
    t = run.time(:);
    if nargin < 3 || isempty(nKnots)
        m = numel(t);
    else
        m = min(nKnots, numel(t));
    end
    if m < 2 || any(diff(t) == 0)
        error('Each run needs at least two distinct time samples.');
    end

    knots = linspace(min(t), max(t), m).';
    x0 = ones(m,1);
    fun = @(x) residual(x, knots, run, dt);
    x = lsqnonlin(fun, x0, 0.85*x0, 1.15*x0, opts);
    s = interp1(knots, x, t, 'linear');

    drag_study_scale('set', knots, x);
    [Rp, Vp, tp] = propagate_traj(run.R(1,:), run.V(1,:), ...
                                t(1), t(end), dt);
    [tp, idx] = sort(tp(:));
    Rp = interp1(tp, Rp(idx,:), t);
    Vp = interp1(tp, Vp(idx,:), t);

    out(j).label = run.label;
    out(j).loftDeg = run.loftDeg;
    out(j).time = t;
    out(j).scale = s;
    out(j).knots = knots;
    out(j).knotScale = x;
    out(j).rmseR = sqrt(mean(sum((Rp-run.R).^2,2)));
    out(j).rmseV = sqrt(mean(sum((Vp-run.V).^2,2)));
end
end

function e = residual(x, knots, run, dt)
t = run.time(:);
drag_study_scale('set', knots, x);
[R, V, tp] = propagate_traj(run.R(1,:), run.V(1,:), ...
                           t(1), t(end), dt);
[tp, idx] = sort(tp(:));
R = interp1(tp, R(idx,:), t);
V = interp1(tp, V(idx,:), t);

n = numel(t);
% Normalize position/velocity errors and penalize irregular scales.
eR = reshape((R-run.R)/100, [], 1)/sqrt(n);
eV = reshape((V-run.V)/1, [], 1)/sqrt(n);
lambda = 1; mu = 0.1;
e = [eR; eV; sqrt(lambda)*diff(x)/0.05/sqrt(numel(x)-1); ...
     sqrt(mu)*(x-1)/0.05/sqrt(numel(x))];
end
