%% Inputs: map these lines to the fields in your workspace truth struct.
r = [truth.x(:),  truth.y(:),  truth.z(:)];     % ECI position, m
v = [truth.vx(:), truth.vy(:), truth.vz(:)];    % ECI velocity, m/s
t = truth.time(:);                              % Truth timestamps, s
% For uniformly sampled truth without timestamps, replace the line above with:
% t = truthStartTime + (0:size(r,1)-1)' * truthSampleInterval;

windows = [30 55 80 105 130 150];  % Or explicit pairs, e.g. [30 55; 80 105]
dt = 0.1;                          % Propagator time step, s
%% Comparison frame (enable one block)
% ECI: use the input positions directly.
toPlotFrame = [];
names = {'ECI X','ECI Y','ECI Z','ECI position magnitude'};

% ACD: comment the two ECI lines above and uncomment both lines below.
% toPlotFrame = @myACD;    % (time, R_ECI) -> N-by-3 ACD in meters
% names = {'Altitude','Crossrange','Downrange','RSS of ACD components'};

%% Generate segments
if isvector(windows)
    bounds = windows(:);
    windows = [bounds(1:end-1),bounds(2:end)];
end
if any(windows(:)<t(1) | windows(:)>t(end))
    error('Window times must be within the truth record [%g, %g] s.',t(1),t(end));
end

seedR = interp1(t,r,windows(:,1),'linear');
seedV = interp1(t,v,windows(:,1),'linear');
segments = repmat(struct('time',[],'R',[],'V',[]),size(windows,1),1);
for k = 1:size(windows,1)
    tStart = windows(k,1);
    tStop  = windows(k,2);
    [rk,vk,tk] = propagate_traj(seedR(k,:),seedV(k,:),tStart,tStop,dt);
    tk = tk(:);
    if isempty(tk)
        error('Segment %d returned no samples for [%g, %g] s.',k,tStart,tStop);
    end
    % Compare scalar endpoints separately. Cap roundoff tolerance well below dt.
    tol = min(dt/1000, max(dt*1e-9,32*eps(max([1;abs(tk);abs(tStart);abs(tStop)]))));
    if abs(tk(1)-tStart)<=tol, tk(1) = tStart; end
    if abs(tk(end)-tStop)<=tol, tk(end) = tStop; end
    % Clip overshoot by interpolating only inside the returned history.
    if tk(1)<tStart || tk(end)>tStop
        times = unique([tk(tk>=tStart & tk<=tStop);tStart;tStop]);
        times = times(times>=tk(1) & times<=tk(end));
        if isempty(times) || numel(tk)<2
            error('Segment %d has insufficient samples in [%g, %g] s.',k,tStart,tStop);
        end
        rk = interp1(tk,rk,times,'linear');
        vk = interp1(tk,vk,times,'linear');
        tk = times;
    end
    if tk(1)~=tStart || tk(end)~=tStop
        warning(['Segment %d requested [%.15g, %.15g] s, covered [%.15g, %.15g] s. ' ...
                 'Start offset %+.6g s, stop shortfall %+.6g s (dt=%g s).'], ...
                 k,tStart,tStop,tk(1),tk(end),tk(1)-tStart,tStop-tk(end),dt);
    end
    segments(k).time = tk;
    segments(k).R = rk;
    segments(k).V = vk;
end

%% Compare truth and propagated positions
truthPos = r;
if ~isempty(toPlotFrame)
    truthPos = toPlotFrame(t,r);
end
errorNames = names;
errorNames{4} = 'RSS position error';
[trajectoryFig,ax] = makeFigure('Trajectory comparison','Position (m)',names);
[errorFig,ex] = makeFigure('Position error','Error (m)',errorNames);
truthColor = [0 .78 .85];
propColor = [.85 0 .65];
lineWidth = 3.5;
truthSeries = [truthPos,vecnorm(truthPos,2,2)];
for j = 1:4
    plot(ax(j),t,truthSeries(:,j),'-','Color',truthColor,'LineWidth',lineWidth,'DisplayName','Truth');
end

edges = unique(windows(:));
for a = [ax ex]
    xline(a,edges,'k--','HandleVisibility','off');
    xlim(a,edges([1 end]));
end
for k = 1:numel(segments)
    ts = segments(k).time;
    propagatedPos = segments(k).R;
    alignedTruthPos = interp1(t,r,ts,'linear');
    if ~isempty(toPlotFrame)
        propagatedPos = toPlotFrame(ts,propagatedPos);
        alignedTruthPos = toPlotFrame(ts,alignedTruthPos);
    end
    delta = propagatedPos-alignedTruthPos;
    y = [propagatedPos,vecnorm(propagatedPos,2,2)];
    e = [delta,vecnorm(delta,2,2)];
    for j = 1:4
        if k>1 && windows(k-1,2)==windows(k,1) && prevTime==ts(1) && ts(1)==windows(k,1)
            plot(ax(j),[prevTime ts(1)],[prevY(j) y(1,j)],'--','Color',propColor, ...
                'LineWidth',lineWidth,'HandleVisibility','off');
            plot(ex(j),[prevTime ts(1)],[prevE(j) e(1,j)],'--','Color',propColor, ...
                'LineWidth',lineWidth,'HandleVisibility','off');
        end
        p = plot(ax(j),ts,y(:,j),'-','Color',propColor, ...
            'LineWidth',lineWidth,'DisplayName','Propagated');
        if k>1, p.HandleVisibility = 'off'; end
        plot(ex(j),ts,e(:,j),'Color',propColor,'LineWidth',lineWidth);
        % Highlight the actual first error sample for each segment.
        plot(ex(j),ts(1),e(1,j),'o','MarkerSize',4,'Color',propColor);
    end
    prevTime = ts(end);
    prevY = y(end,:);
    prevE = e(end,:);
end
legend(ax(1),'Location','best');
linkaxes(ax,'x');
linkaxes(ex,'x');

function [f,ax] = makeFigure(name,units,names)
f = figure('Name',name,'Color','w');
tl = tiledlayout(f,3,2,'TileSpacing','compact','Padding','compact');
ax = gobjects(1,4);
for j = 1:3
    ax(j) = nexttile(tl,2*j-1);
end
ax(4) = nexttile(tl,2,[3 1]);
for j = 1:4
    set(ax(j),'FontSize',13,'FontWeight','bold','LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1);
    hold(ax(j),'on');
    grid(ax(j),'on');
    title(ax(j),names{j});
    ylabel(ax(j),units);
end
xlabel(ax(3),'Time (s)');
xlabel(ax(4),'Time (s)');
title(tl,name);
end
