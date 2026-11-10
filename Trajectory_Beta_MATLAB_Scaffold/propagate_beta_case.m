function sim = propagate_beta_case(R_0, V_0, t_0, t_final, dt, drag)
%PROPAGATE_BETA_CASE Adapter for the production trajectory propagator.
%
% TODO: connect this adapter to your production function. One possible
% extension of the existing interface is:
% [R,V,time,acc_grav,acc_drag] = ...
%     propagate_traj(R_0,V_0,t_0,t_final,dt,drag);
%
% Use the existing predictor's position-based gravity helper. In the drag
% helper, use the current altitude calculation, Atmos, and Earth rotation
% to obtain the air-relative Mach number. Replace the hard-coded lookup:
% beta = interp1(Mtable, coefficientTable, Mach);
% with the configured lookup at that force evaluation's stage time t:
% query = struct('Mach', Mach);
% beta = drag.betaAt(t, query, segmentMode);
% Supply any additional verified table coordinates in query. Thread t and
% the drag context through the caller if the helper currently takes R and V.
% For a point query use [] as segmentMode. Split crossing steps at times in
% drag.switchTimes and hold the selected mode through each segment's stages.
%
% TODO: confirm the units and normalization of beta in the existing drag
% acceleration expression and use the lookup output in that expression.
% TODO: map the production outputs to the README history contract. After
% unit, frame, and shape mapping, the basic output assignments are:
% sim.time = time(:);
% sim.signals = struct('R', R, 'V', V, 'acc_grav', acc_grav, 'acc_drag', acc_drag);

error('BetaScaffold:TODO', ...
    'Connect propagate_beta_case.m to your production propagator. See README.md.');
end
