function aL = lift_accel(R, V, t, aero)
%LIFT_ACCEL Calibrated lift acceleration, perpendicular to relative wind.
% Positive CL acts toward local radial-up, projected perpendicular to wind.
% This is a prescribed lift-up direction, NOT a body-attitude/roll model.
% Reuses dynamic pressure from Drag; does not use table-derived CL.

aL = zeros(size(V));
if isempty(aero)
    return
end

persistent K_norm
if isempty(K_norm)
    K_norm = evalin('base','K_norm'); % Same normalization as Drag.m
end

R = reshape(R,1,3);
V = reshape(V,1,3);
omegaEarth = [0 0 7e-5];       % Match cfg.omegaEarth in Drag.m
Vrel = V - cross(omegaEarth,R);
vhat = Vrel/norm(Vrel);
rhat = R/norm(R);

liftDir = rhat - dot(rhat,vhat)*vhat;
L = norm(liftDir);
if L < 1e-8                  % Undefined lift-up direction for radial flight
    return
end

cL = lift_study_coefficient(t);
aL = reshape(aero.qbar * K_norm * cL * liftDir/L, size(aL));
end
