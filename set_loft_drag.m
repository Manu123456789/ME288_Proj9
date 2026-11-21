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

How to :
Exactly. Now that you have the optimized drag-scale histories for each loft angle, the next step is to build a **two-dimensional interpolation** over loft angle and time.

You do **not** need to run `lsqnonlin` again. We're simply interpolating the previously optimized results.

For example, with fitted trajectories at 15° and 20°, the scale at 16.7° would be

\\[ \boxed{s(16.7^\circ,\tau)=0.66s(15^\circ,\tau)+0.34s(20^\circ,\tau)} \\]

where \\(\tau\\) is normalized flight time. This accounts for the possibility that the trajectories have different durations.

## **1. New MATLAB function**

Download set_loft_drag.m

This uses your existing `out(j).loftDeg`, `out(j).knots`, and `out(j).knotScale`, and automatically loads the interpolated correction into `drag_study_scale`.

```

function [tNodes, scale] = set_loft_drag(out, loftDeg, t0, tf)
% Interpolate drag correction between fitted loft angles.
% out = output from fit_loft_drag_time

[lofts, order] = sort([out.loftDeg]);

if any(diff(lofts) <= 0) || ...
        loftDeg < lofts(1) || loftDeg > lofts(end)
    error('Loft outside fitted range or duplicate lofts.');
end

if t0 == tf
    error('Start and end times must differ.');
end

% Find neighboring lofts
i = find(lofts <= loftDeg, 1, 'last');
j = find(lofts >= loftDeg, 1, 'first');

a = out(order(i));
b = out(order(j));

% Normalize time for different trajectory durations
tauA = (a.knots(:)-min(a.time)) / ...
       (max(a.time)-min(a.time));
tauB = (b.knots(:)-min(b.time)) / ...
       (max(b.time)-min(b.time));

% Combine optimization time nodes
tau = unique([tauA; tauB]);

scaleA = interp1(tauA,a.knotScale(:),tau);
scaleB = interp1(tauB,b.knotScale(:),tau);

% Interpolate between loft angles
if i == j
    scale = scaleA;
else
    w = (loftDeg-lofts(i))/(lofts(j)-lofts(i));
    scale = (1-w)*scaleA + w*scaleB;
end

% Convert to the new propagation time range
tNodes = t0 + (tf-t0)*tau;

% Apply correction through existing Drag.m
drag_study_scale('set',tNodes,scale);

end

```

## **2. How to use it in your propagator**

Your main script becomes very simple.

```
% New loft angle
loft = 16.7;

% Desired propagation time range
t0 = 0;
tf = 120;  % Example: 120 seconds

% Load interpolated drag correction
[tScale, scale] = set_loft_drag(out,loft,t0,tf);

% Run existing propagator
[R,V,time] = propagate_traj(R0,V0,t0,tf,dt);

% Retrieve correction at every propagated time step
sHistory = arrayfun(@drag_study_scale,time);

% Reset to nominal drag
drag_study_scale('reset');

% Plot interpolated correction
figure
plot(time,sHistory)
xlabel('Time (s)')
ylabel('Drag multiplier')
grid on
```

Your existing `Drag.m` continues using:

```
inverseBeta = drag_study_scale(t) * model.K_norm * CD;
```

**No modifications are needed to the propagator, drag function, or scale helper**, provided they already support the time-dependent corrections from our previous implementation.

## **3. What this gives you**

Illustrative interpolated correction histories

Illustrative fitted 15° and 20° correction profiles and the interpolated profile at 16.7°.

15° fitted scale16.7° interpolated scale20° fitted scale

0.940.9650.991.0151.040%20%40%60%80%100%

Example only — these are not your fitted coefficients.

You now have a lookup procedure that accepts any loft angle within your fitted range and returns a time-varying drag multiplier. Different calibration trajectories can have different knot locations and durations.

**One important assumption:** Normalizing time means we're matching the same *fraction of the flight*, not necessarily the same physical time or Mach/Reynolds conditions. That is a reasonable first approach, but it must be validated on loft angles that were not used in fitting. If flight durations or aerodynamic regimes differ substantially, this interpolation may not generalize well.
