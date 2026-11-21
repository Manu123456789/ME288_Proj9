% Enter the truth trajectories in propagation order; R,V are Nx3 ECI.
makeRun = @(label,loft,t,R,V) struct('label',label, 'loftDeg',loft, ...
                                     'time',t(:), 'R',R, 'V',V);
runs(1) = makeRun('Case A', 10, tA, RA, VA);
runs(2) = makeRun('Case B', 15, tB, RB, VB);

% One variable per truth time (potentially expensive):
out = fit_loft_drag_time(runs, dt);
% For faster fitting, optimize 25 nodes, return scale at EVERY saved time:
% out = fit_loft_drag_time(runs, dt, 25);

figure
hold on
for j = 1:numel(out)
    plot(out(j).time, out(j).scale, 'DisplayName',out(j).label);
end
xlabel('Time (s)'); ylabel('Drag multiplier'); legend('show'); grid on

% Example: replay the optimized schedule for run 1
k = 1;
drag_study_scale('set', out(k).knots, out(k).knotScale);
[Rfit,Vfit,tfit] = propagate_traj(runs(k).R(1,:),runs(k).V(1,:), ...
                                  runs(k).time(1),runs(k).time(end),dt);
drag_study_scale('reset');
