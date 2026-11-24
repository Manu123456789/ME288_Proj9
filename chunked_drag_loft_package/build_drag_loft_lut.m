function lut = build_drag_loft_lut(out,runs,dt,source)
% Convert time-calibrated drag to f(loft,Mach,Re).
% Default uses truth states: build_drag_loft_lut(out,runs)
% Optional replay: build_drag_loft_lut(out,runs,dt,'propagated')
% Keep Drag.m in time-based mode during a propagated replay.
if nargin < 4 || isempty(source), source = 'truth'; end
if ~ismember(source,{'truth','propagated'})
    error('source must be truth or propagated.');
end
if numel(out) ~= numel(runs) || isempty(out)
    error('out and runs must contain the same nonempty cases.');
end
if strcmp(source,'propagated') && (nargin < 3 || isempty(dt))
    error('Provide the original integration step/output spacing for replay.');
end
cleanup = onCleanup(@() drag_study_scale('reset'));
drag_study_scale('reset');
drag_loft_state_scale('reset');
blocks = cell(numel(out),1);
for j = 1:numel(out)
    a = out(j); r = runs(j);
    if a.loftDeg ~= r.loftDeg || ...
            (isfield(a,'label') && isfield(r,'label') && ~strcmp(a.label,r.label))
        error('Case %d does not match between out and runs.',j);
    end
    t = a.time(:);
    if numel(t)<2 || ~(all(diff(t)>0) || all(diff(t)<0))
        error('Case %d requires strictly monotonic fit times.',j);
    end
    [tk,ix] = sort(a.knots(:));
    sk = a.knotScale(:);
    if numel(sk) ~= numel(tk), error('Case %d has mismatched knots and scales.',j); end
    sk = sk(ix);
    if numel(tk)<2 || any(diff(tk)<=0) || any(~isfinite(sk)) || any(sk<=0)
        error('Case %d has invalid knots or nonpositive drag scales.',j);
    end
    s = interp1(tk,sk,t,'linear',NaN);
    [tr,ix] = sort(r.time(:));
    if any(diff(tr)<=0) || numel(tr)<2
        error('Case %d has invalid truth timestamps.',j);
    end
    R = interp1(tr,r.R(ix,:),t,'linear',NaN);
    V = interp1(tr,r.V(ix,:),t,'linear',NaN);
    if any(~isfinite([R(:);V(:);s]))
        error('Case %d has missing data or fit times outside truth/knots.',j);
    end
    if strcmp(source,'propagated')
        drag_study_scale('set',tk,sk);
        [Rp,Vp,tp] = propagate_traj(R(1,:),V(1,:),t(1),t(end),dt);
        [tp,ix] = sort(tp(:));
        R = interp1(tp,Rp(ix,:),t,'linear',NaN);
        V = interp1(tp,Vp(ix,:),t,'linear',NaN);
        if any(~isfinite([R(:);V(:)]))
            error('Replay %d does not cover all fitted timestamps.',j);
        end
    end
    M = nan(size(t)); Re = M;
    for k = 1:numel(t)
        [~,aero] = Drag(R(k,:),V(k,:),t(k));
        if ~isempty(aero)
            M(k) = aero.Mach;
            Re(k) = aero.Re;
        end
    end
    keep = isfinite(M) & isfinite(Re) & M>=0 & Re>=0;
    if ~any(keep), error('Case %d has no valid aero samples.',j); end
    if any(~keep), warning('Skipped %d invalid aero samples in case %d.',sum(~keep),j); end
    n = sum(keep);
    blocks{j} = table(repmat(j,n,1),t(keep),repmat(a.loftDeg,n,1), ...
        M(keep),Re(keep),s(keep), 'VariableNames', ...
        {'Case','Time','Loft','Mach','Re','Scale'});
    drag_study_scale('reset');
end
lut = make_drag_loft_lut(vertcat(blocks{:}));
lut.source = source;
end
