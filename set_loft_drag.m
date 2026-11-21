function [tNodes, scale] = set_loft_drag(out, loftDeg, t0, tf)
%SET_LOFT_DRAG Interpolate fitted time-varying drag scale between lofts.
% out is the struct array returned by fit_loft_drag_time.
% t0 and tf are the start and end times for the new propagation.

[lofts, order] = sort([out.loftDeg]);
if any(diff(lofts) <= 0) || loftDeg < lofts(1) || loftDeg > lofts(end)
    error('Use one fit per loft and a loft within the fitted range.');
end
if t0 == tf
    error('Start and end times must differ.');
end

i = find(lofts <= loftDeg, 1, 'last');
j = find(lofts >= loftDeg, 1, 'first');
a = out(order(i));
b = out(order(j));

% Normalize each fitted time axis to [0,1] to handle different durations.
tauA = (a.knots(:) - min(a.time)) / (max(a.time) - min(a.time));
tauB = (b.knots(:) - min(b.time)) / (max(b.time) - min(b.time));
tau = unique([tauA; tauB]);

scaleA = interp1(tauA, a.knotScale(:), tau, 'linear');
scaleB = interp1(tauB, b.knotScale(:), tau, 'linear');

if i == j
    scale = scaleA;
else
    w = (loftDeg - lofts(i)) / (lofts(j) - lofts(i));
    scale = (1 - w)*scaleA + w*scaleB;
end

% Install the time-varying scale for the existing Drag(R,V,t).
tNodes = t0 + (tf - t0)*tau;
drag_study_scale('set', tNodes, scale);
end
