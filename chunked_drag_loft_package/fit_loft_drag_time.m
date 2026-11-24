function out = fit_loft_drag_time(runs, dt, nKnots)
% Fit time-varying drag scales for each truth trajectory.
% runs(j): label, loftDeg, time (Nx1), R (Nx3), V (Nx3)

opts = optimoptions('lsqnonlin','Display','off');
cleanup = onCleanup(@() drag_study_scale('reset'));
out = struct([]);

% >>> ADDED: adaptive-refinement settings (change these as needed)
maxKnots = 60;
refinePasses = 2;
newPerPass = 15;
% <<< ADDED

for j = 1:numel(runs)
    run = runs(j);
    t = run.time(:);

    if nargin < 3 || isempty(nKnots)
        % >>> CHANGED: start with a coarse fit, not one knot per sample
        m = min(20,numel(t));
        % <<< CHANGED
    else
        m = min(nKnots,numel(t));
    end

    if m < 2 || any(diff(t) == 0) || ...
            ~(all(diff(t) > 0) || all(diff(t) < 0))
        error('Time samples must be strictly monotonic and contain at least two points.');
    end

    knots = linspace(min(t), max(t), m).';
    x0 = ones(m,1);

    fun = @(x) residual(x, knots, run, dt);
    x = lsqnonlin(fun, x0, 0.85*x0, 1.15*x0, opts);

    % >>> ADDED: refine knots where velocity error magnitude changes fastest
    for pass = 1:refinePasses
        if numel(knots) >= maxKnots
            break
        end
        drag_study_scale('set', knots, x);
        [~,Vp,tp] = propagate_traj( ...
            run.R(1,:), run.V(1,:), t(1), t(end), dt);
        [tp,idx] = sort(tp(:));
        Vp = interp1(tp,Vp(idx,:),t);
        vErr = sqrt(sum((Vp-run.V).^2,2));
        growth = max(diff(vErr)./abs(diff(t)),0);
        growth = movmean(growth,5);
        [~,order] = sort(growth,'descend');

        oldKnots = knots;
        oldX = x;
        minSpacing = 0.25*min(diff(oldKnots));
        count = 0;
        for q = order(:).'
            candidate = (t(q)+t(q+1))/2;
            if min(abs(knots-candidate)) > minSpacing
                knots = sort([knots; candidate]);
                count = count+1;
            end
            if count >= min(newPerPass, maxKnots-numel(oldKnots))
                break
            end
        end
        if count == 0
            break
        end

        x0 = interp1(oldKnots,oldX,knots,'linear'); % warm start
        fun = @(z) residual(z, knots, run, dt);
        x = lsqnonlin(fun,x0,0.85*ones(size(x0)), ...
                     1.15*ones(size(x0)),opts);
    end
    % <<< ADDED

    s = interp1(knots, x, t, 'linear');

    drag_study_scale('set', knots, x);
    [Rp,Vp,tp] = propagate_traj( ...
        run.R(1,:), run.V(1,:), t(1), t(end), dt);
    [tp,idx] = sort(tp(:));
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
[R,V,tp] = propagate_traj( ...
    run.R(1,:),run.V(1,:),t(1),t(end),dt);
[tp,idx] = sort(tp(:));
R = interp1(tp,R(idx,:),t);
V = interp1(tp,V(idx,:),t);

n = numel(t);
eR = reshape((R-run.R)/100,[],1)/sqrt(n);
eV = reshape((V-run.V)/1,[],1)/sqrt(n);

lambda = 1;
mu = 0.1;
e = [eR; eV; ...
     sqrt(lambda)*diff(x)/0.05/sqrt(numel(x)-1); ...
     sqrt(mu)*(x-1)/0.05/sqrt(numel(x))];
end
