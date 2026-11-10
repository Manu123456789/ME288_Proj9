%% BASELINE COMPARISON: existing propagator vs. truth trajectory
% No beta-model changes. Requires propagate_traj.m on the MATLAB path.

%% INPUTS (edit only this section)
truthFile = 'truth_data.mat';  % MAT file containing the truth struct
truthName = 'truth';           % Name of the struct stored in that file

fieldTime = 'time';            % Truth time field: N x 1
fieldR    = 'R';               % Truth position field: N x 3
fieldV    = 'V';               % Truth velocity field: N x 3

t1 = 0;                       % Initial time [s]
t2 = -20;                     % Final time [s]; may be before t1
dt = 0.01;                    % Propagator step magnitude [s]

%% LOAD AND CHECK TRUTH
S = load(truthFile);
truth = S.(truthName);
tTruth = truth.(fieldTime); tTruth = tTruth(:);
RTruth = truth.(fieldR);
VTruth = truth.(fieldV);

assert(size(RTruth,1)==numel(tTruth) && size(RTruth,2)==3, ...
    'Truth position must be N-by-3, with N truth timestamps.');
assert(size(VTruth,1)==numel(tTruth) && size(VTruth,2)==3, ...
    'Truth velocity must be N-by-3, with N truth timestamps.');
assert(all(isfinite(tTruth)) && all(isfinite(RTruth(:))) && ...
    all(isfinite(VTruth(:))), 'Truth contains NaN or Inf.');
assert(isscalar(dt) && isfinite(dt) && dt>0 && t1~=t2, ...
    'Use dt > 0 and different finite start/end times.');

% Sort truth once, so interpolation works even if the file is reversed.
[tTruth, order] = sort(tTruth);
RTruth = RTruth(order,:);
VTruth = VTruth(order,:);
assert(all(diff(tTruth)>0), 'Truth timestamps must be unique.');
assert(min(t1,t2)>=tTruth(1) && max(t1,t2)<=tTruth(end), ...
    'Requested interval lies outside the truth history.');

%% SEED AND RUN EXISTING PROPAGATOR
R0 = interp1(tTruth, RTruth, t1, 'linear');
V0 = interp1(tTruth, VTruth, t1, 'linear');
[R, V, tModel, ~, ~] = propagate_traj(R0, V0, t1, t2, dt);
tModel = tModel(:);

assert(size(R,1)==numel(tModel) && size(R,2)==3 && ...
       size(V,1)==numel(tModel) && size(V,2)==3, ...
       'Propagator outputs R and V must be N-by-3.');
tol = 1e-8 * max([1,abs(t1),abs(t2)]);
assert(abs(tModel(1)-t1)<tol && abs(tModel(end)-t2)<tol, ...
    'Propagator did not cover the requested interval exactly.');

%% COMPUTE ERRORS AT THE PROPAGATOR TIME STAMPS
RRef = interp1(tTruth, RTruth, tModel, 'linear');
VRef = interp1(tTruth, VTruth, tModel, 'linear');
errors = {R-RRef, V-VRef};  % Signed component errors, model minus truth
names  = {'Position', 'Velocity'};
units  = {'m', 'm/s'};
axes3  = {'X','Y','Z'};

%% TWO FIGURES: 3 COMPONENT PLOTS LEFT, TOTAL MAGNITUDE RIGHT
for k = 1:2
    E = errors{k};
    totalError = sqrt(sum(E.^2,2));  % 3D Euclidean norm per time sample

    figure('Name',[names{k} ' baseline error'],'Color','w');
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
    for j = 1:3
        nexttile(2*j-1);  % Left column: tiles 1, 3, 5
        plot(tModel,E(:,j),'LineWidth',1.2);
        ylabel(sprintf('%s error [%s]',axes3{j},units{k}));
        grid on;
        if j==3, xlabel('Time [s]'); end
    end

    nexttile(2,[3 1]);  % Right column: one plot spanning all three rows
    plot(tModel,totalError,'LineWidth',1.5); hold on;
    plot(tModel(end),totalError(end),'o','MarkerSize',7,'LineWidth',1.5);
    grid on;
    xlabel('Time [s]');
    ylabel(sprintf('3D error magnitude [%s]',units{k}));
    title('Total error (endpoint marked)');
    sgtitle([names{k} ' error: propagator minus truth']);

    fprintf('\n%s error [%s]\n',names{k},units{k});
    fprintf('  Final XYZ:   [% .6g, % .6g, % .6g]\n',E(end,:));
    fprintf('  Final total: %.6g\n',totalError(end));
    fprintf('  RMS total:   %.6g\n',sqrt(mean(totalError.^2)));
    fprintf('  Peak total:  %.6g\n',max(totalError));
end
