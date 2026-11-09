function sim = propagate_beta_case(R_0, V_0, t_0, t_final, dt, drag)
%PROPAGATE_BETA_CASE Adapter for the production trajectory propagator.
%
% TODO: connect this adapter to your production function. One possible
% extension of the existing interface is:
% [R,V,time,acc_grav,acc_drag] = ...
%     propagate_traj(R_0,V_0,t_0,t_final,dt,drag);
%
% At each force evaluation, compute query using existing flow calculations.
% Evaluate beta with drag.betaAt(t, query, segmentMode). Use [] for a point
% query. At each time in drag.switchTimes, end the crossing step and keep
% segmentMode fixed for the stages on each side of the event.
%
% TODO: confirm how beta enters the existing drag acceleration equation.
% TODO: map the production outputs to the README history contract. After
% unit, frame, and shape mapping, the basic output assignments are:
% sim.time = time(:);
% sim.signals = struct('R', R, 'V', V, 'acc_grav', acc_grav, 'acc_drag', acc_drag);

error('BetaScaffold:TODO', ...
    'Connect propagate_beta_case.m to your production propagator. See README.md.');
end
