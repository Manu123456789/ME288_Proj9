% Requires your existing runs, dt, Drag.m and propagate_traj.
% BUILD: keep the original time-based multiplier in Drag.m active.
chunkSeconds = 30;
nKnots = 8;   % initial knots per chunk; original adaptive settings still apply
[out,chunkRuns] = fit_loft_drag_chunks(runs,dt,chunkSeconds,nKnots);

% Use the corresponding CHUNK histories, not the unsplit runs.
lut = build_drag_loft_lut(out,chunkRuns);
% Optional: associate scales with the corrected propagated states instead.
% lut = build_drag_loft_lut(out,chunkRuns,dt,'propagated');
save('drag_loft_chunked.mat','lut','out','chunkRuns');

figure
scatter3(lut.samples.Mach,lut.samples.Re,lut.samples.Loft, ...
    12,lut.samples.Scale,'filled');
xlabel('Mach'); ylabel('Re (table scale)'); zlabel('Loft (deg)');
colorbar; grid on; title('Chunk-based calibration samples');

% RUNTIME: switch Drag.m to its state-based multiplier (see README).
% Run the following separately after that change.
% load('drag_loft_chunked.mat','lut');
% drag_loft_state_scale('set',lut,16.7);
% [R,V,time] = propagate_traj(R0,V0,t0,tf,dt);
% drag_loft_state_scale('reset');
