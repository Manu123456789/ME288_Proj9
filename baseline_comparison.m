
%% BASELINE COMPARISON: Existing Propagator vs. Truth Trajectory
% Requires your existing propagate_traj.m on the MATLAB path.

%% =========================================================
%  INPUTS (EDIT THIS SECTION)
% ==========================================================

truthFile = 'truth_data.mat';   % MAT file containing truth trajectory struct
truthName = 'truth';            % Struct variable name in MAT file

fieldTime = 'time';             % Truth time field (N x 1)
fieldR    = 'R';                % Truth position field (N x 3)
fieldV    = 'V';                % Truth velocity field (N x 3)

t1 = 0;                        % Initial time [s]
t2 = -20;                      % Final time [s] (forward or backward)
dt = 0.01;                     % Propagator time step [s]


%% =========================================================
%  LOAD TRUTH AND EXTRACT INITIAL STATE
% ==========================================================

S = load(truthFile);
truth = S.(truthName);

tTruth = truth.(fieldTime);
tTruth = tTruth(:);

RTruth = truth.(fieldR);
VTruth = truth.(fieldV);

% Reverse a descending truth history for interpolation
if tTruth(1) > tTruth(end)
    tTruth = flipud(tTruth);
    RTruth = flipud(RTruth);
    VTruth = flipud(VTruth);
end

% Essential checks
assert(all(diff(tTruth) > 0), ...
    'Truth timestamps must be ordered and unique.');

assert(min(t1,t2) >= tTruth(1) && max(t1,t2) <= tTruth(end), ...
    'Requested time interval must be within the truth history.');

% Extract initial position and velocity at t1
R0 = interp1(tTruth, RTruth, t1, 'linear');
V0 = interp1(tTruth, VTruth, t1, 'linear');


%% =========================================================
%  RUN EXISTING PROPAGATOR
% ==========================================================

[R, V, tModel, ~, ~] = ...
    propagate_traj(R0, V0, t1, t2, dt);

tModel = tModel(:);

tol = 1e-8 * max([1, abs(t1), abs(t2)]);

assert(abs(tModel(1)-t1) < tol && abs(tModel(end)-t2) < tol, ...
    'Propagator did not reach the requested start and end times.');

% Convert model history to increasing time for interpolation
if tModel(1) > tModel(end)
    tModelUp = flipud(tModel);
    RUp = flipud(R);
    VUp = flipud(V);
else
    tModelUp = tModel;
    RUp = R;
    VUp = V;
end

% Remove endpoint roundoff
tModelUp([1 end]) = [min(t1,t2), max(t1,t2)];


%% =========================================================
%  ALIGN TRUTH AND MODEL TIME HISTORIES
% ==========================================================

% Use the coarser time grid for comparison.
% Interpolate the other history to the comparison times.

if median(diff(tTruth)) >= median(diff(tModelUp))

    % Truth has coarser sampling
    tCompare = tTruth(tTruth >= min(t1,t2) & ...
                      tTruth <= max(t1,t2));

else

    % Propagator has coarser sampling
    tCompare = tModelUp;

end

% Include exact initial and final times
tCompare = unique([tCompare; t1; t2]);

% Preserve propagation direction
if t2 < t1
    tCompare = flipud(tCompare);
end

% Interpolate truth onto comparison times
RRef = interp1(tTruth, RTruth, tCompare, 'linear');
VRef = interp1(tTruth, VTruth, tCompare, 'linear');

% Interpolate propagator onto comparison times
RSim = interp1(tModelUp, RUp, tCompare, 'linear');
VSim = interp1(tModelUp, VUp, tCompare, 'linear');


%% =========================================================
%  COMPUTE ERRORS
% ==========================================================

% Signed component errors: propagated - truth
dR = RSim - RRef;
dV = VSim - VRef;

% Total 3D error magnitudes
posError = vecnorm(dR, 2, 2);
velError = vecnorm(dV, 2, 2);

% Final errors
finalPosError = posError(end);
finalVelError = velError(end);


%% =========================================================
%  FIGURE 1: POSITION ERRORS
% ==========================================================

figure('Name','Position Baseline Error','Color','w');

tiledlayout(3,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

xyz = {'X','Y','Z'};

% Left column: X, Y, Z component errors
for j = 1:3

    nexttile(2*j-1);

    plot(tCompare, dR(:,j), 'LineWidth',1.2);

    ylabel(sprintf('%s Error [m]',xyz{j}));
    grid on;

    if j == 3
        xlabel('Time [s]');
    end

end

% Right column: total position error magnitude
nexttile(2,[3 1]);

plot(tCompare, posError, 'LineWidth',1.5);
hold on;

plot(tCompare(end), finalPosError, 'o', ...
    'MarkerSize',7, 'LineWidth',1.5);

xlabel('Time [s]');
ylabel('3D Position Error [m]');

title(sprintf('Final Position Error: %.6g m',finalPosError));

grid on;

sgtitle('Position Error: Propagator - Truth');


%% =========================================================
%  FIGURE 2: VELOCITY ERRORS
% ==========================================================

figure('Name','Velocity Baseline Error','Color','w');

tiledlayout(3,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% Left column: X, Y, Z component errors
for j = 1:3

    nexttile(2*j-1);

    plot(tCompare, dV(:,j), 'LineWidth',1.2);

    ylabel(sprintf('%s Error [m/s]',xyz{j}));
    grid on;

    if j == 3
        xlabel('Time [s]');
    end

end

% Right column: total velocity error magnitude
nexttile(2,[3 1]);

plot(tCompare, velError, 'LineWidth',1.5);
hold on;

plot(tCompare(end), finalVelError, 'o', ...
    'MarkerSize',7, 'LineWidth',1.5);

xlabel('Time [s]');
ylabel('3D Velocity Error [m/s]');

title(sprintf('Final Velocity Error: %.6g m/s',finalVelError));

grid on;

sgtitle('Velocity Error: Propagator - Truth');


%% =========================================================
%  PRINT ERROR STATISTICS
% ==========================================================

fprintf('\n========================================\n');
fprintf('       BASELINE PROPAGATION ERRORS\n');
fprintf('========================================\n');

fprintf('Initial Time: %.6f s\n', t1);
fprintf('Final Time:   %.6f s\n', t2);
fprintf('Time Step:    %.6f s\n', dt);

% Position
fprintf('\n--- POSITION ERRORS [m] ---\n');

fprintf('Final X Error: % .6g\n', dR(end,1));
fprintf('Final Y Error: % .6g\n', dR(end,2));
fprintf('Final Z Error: % .6g\n', dR(end,3));

fprintf('FINAL 3D POSITION ERROR: %.6g m\n', finalPosError);
fprintf('RMS 3D Position Error:   %.6g m\n', sqrt(mean(posError.^2)));
fprintf('Peak 3D Position Error:  %.6g m\n', max(posError));

% Velocity
fprintf('\n--- VELOCITY ERRORS [m/s] ---\n');

fprintf('Final X Error: % .6g\n', dV(end,1));
fprintf('Final Y Error: % .6g\n', dV(end,2));
fprintf('Final Z Error: % .6g\n', dV(end,3));

fprintf('FINAL 3D VELOCITY ERROR: %.6g m/s\n', finalVelError);
fprintf('RMS 3D Velocity Error:   %.6g m/s\n', sqrt(mean(velError.^2)));
fprintf('Peak 3D Velocity Error:  %.6g m/s\n', max(velError));

% Endpoint comparison
fprintf('\n--- FINAL POSITION COMPARISON [m] ---\n');

fprintf('Truth:      [% .6g, % .6g, % .6g]\n', RRef(end,:));
fprintf('Propagated: [% .6g, % .6g, % .6g]\n', RSim(end,:));

fprintf('\n========================================\n');
