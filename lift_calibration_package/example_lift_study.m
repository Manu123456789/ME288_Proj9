% Existing inputs:
% runs(j).label, .loftDeg, .time, .R, .V
% outDrag = your previous fit_loft_drag_time(runs,dt) result.
% Run with your existing Drag.m, Atmos.m and propagate_traj.m on the path.

outLift = fit_loft_lift_time(runs,outDrag,dt,12);

figure; hold on
for j = 1:numel(outLift)
    plot(outLift(j).time,outLift(j).CL,'DisplayName',outLift(j).label)
end
xlabel('Time (s)'); ylabel('Fitted effective C_L');
legend('show'); grid on

% Replay an existing fitted loft with fixed drag and fitted lift.
k = 1;
drag_study_scale('set',outDrag(k).knots,outDrag(k).knotScale);
lift_study_coefficient('set',outLift(k).knots,outLift(k).knotCL);
[R,V,t] = propagate_traj(runs(k).R(1,:),runs(k).V(1,:), ...
    runs(k).time(1),runs(k).time(end),dt);
drag_study_scale('reset');
lift_study_coefficient('reset');

% For an interpolated loft, use BOTH front-end setters before propagation:
% set_loft_drag(outDrag,16.7,t0,tf);
% set_loft_lift(outLift,16.7,t0,tf);
% [R,V,t] = propagate_traj(R0,V0,t0,tf,dt);
% drag_study_scale('reset'); lift_study_coefficient('reset');
