function out = fit_loft_lift_time(runs, dragOut, dt, nKnots)
%FIT_LOFT_LIFT_TIME Fit time-varying effective lift CL for each loft.
% runs and dragOut must correspond in order (same cases).
% Drag multipliers stay FIXED while lift is fitted.

if nargin < 4 || isempty(nKnots), nKnots = 12; end
opts = optimoptions('lsqnonlin','Display','iter', ...
    'FiniteDifferenceStepSize',1e-3,'MaxIterations',150);
cleanup = onCleanup(@resetModels); %#ok<NASGU>
out = struct([]);

for j = 1:numel(runs)
    run = runs(j);
    t = run.time(:);
    if numel(t) < 2 || ~(all(diff(t) > 0) || all(diff(t) < 0))
        error('Each run needs at least two monotonic timestamps.');
    end
    if run.loftDeg ~= dragOut(j).loftDeg
        error('runs and dragOut must have matching lofts in matching order.');
    end

    % Reuse already-calibrated drag; optimize ONLY lift.
    drag_study_scale('set',dragOut(j).knots,dragOut(j).knotScale);
    m = max(2,min(nKnots,numel(t)));

    % Baseline: fitted drag only, no lift.
    lift_study_coefficient('reset');
    [Rbase,Vbase,tbase] = propagate_traj( ...
        run.R(1,:),run.V(1,:),t(1),t(end),dt);
    [tbase,idxBase] = sort(tbase(:));
    Rbase = interp1(tbase,Rbase(idxBase,:),t);
    Vbase = interp1(tbase,Vbase(idxBase,:),t);

    knots = linspace(min(t),max(t),m).';
    x0 = zeros(m,1);       % CL = 0 gives your current drag-only model
    lb = -0.2*ones(m,1);  % Example bounds: allow both lift directions
    ub =  0.2*ones(m,1);

    fun = @(x) residual(x,knots,run,dt);
    [x,J,~,exitflag,info] = lsqnonlin(fun,x0,lb,ub,opts);

    lift_study_coefficient('set',knots,x);
    [Rp,Vp,tp] = propagate_traj( ...
        run.R(1,:),run.V(1,:),t(1),t(end),dt);
    [tp,idx] = sort(tp(:));
    Rp = interp1(tp,Rp(idx,:),t);
    Vp = interp1(tp,Vp(idx,:),t);

    out(j).label = run.label;
    out(j).loftDeg = run.loftDeg;
    out(j).time = t;
    out(j).CL = interp1(knots,x,t,'linear');
    out(j).knots = knots;
    out(j).knotCL = x;
    out(j).rmseR0 = sqrt(mean(sum((Rbase-run.R).^2,2)));
    out(j).rmseV0 = sqrt(mean(sum((Vbase-run.V).^2,2)));
    out(j).rmseR = sqrt(mean(sum((Rp-run.R).^2,2)));
    out(j).rmseV = sqrt(mean(sum((Vp-run.V).^2,2)));
    out(j).resnorm = J;
    out(j).exitflag = exitflag;
    out(j).message = info.message;
end
end

function e = residual(x,knots,run,dt)
t = run.time(:);
lift_study_coefficient('set',knots,x);

[R,V,tp] = propagate_traj( ...
    run.R(1,:),run.V(1,:),t(1),t(end),dt);
[tp,idx] = sort(tp(:));
R = interp1(tp,R(idx,:),t);
V = interp1(tp,V(idx,:),t);

n = numel(t);
eR = reshape((R-run.R)/100,[],1)/sqrt(n);
eV = reshape((V-run.V)/1,[],1)/sqrt(n);
lambda = 0.1;  % Penalize rapid changes in CL
mu = 0.01;     % Keep effective CL near zero

e = [eR; eV; ...
     sqrt(lambda)*diff(x)/0.05/sqrt(numel(x)-1); ...
     sqrt(mu)*x/0.05/sqrt(numel(x))];
end

function resetModels()
lift_study_coefficient('reset');
drag_study_scale('reset');
end
