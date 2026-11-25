%% Inputs: map these lines to the fields in your workspace truth struct.
r = [truth.x(:),  truth.y(:),  truth.z(:)];     % ECI position, m
v = [truth.vx(:), truth.vy(:), truth.vz(:)];    % ECI velocity, m/s
t = truth.time(:);                              % Truth timestamps, s
% For uniformly sampled truth without timestamps, replace the line above with:
% t = truthStartTime + (0:size(r,1)-1)' * truthSampleInterval;

windows = [30 55 80 105 130 150];  % Or explicit pairs, e.g. [30 55; 80 105]
dt = 0.1;                          % Propagator time step, s
convert = @myACD;                  % (time, R_ECI) -> [altitude,crossrange,downrange], m

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
    % Snap endpoint roundoff only. Keep actual timestamps for shorter outputs.
    requested = [tStart;tStop];
    if numel(tk)>=2 && all(abs(tk([1 end])-requested)<=16*eps(max(1,abs(requested))))
        tk([1 end]) = requested;
    else
        warning('Segment %d did not return both requested endpoints [%g, %g] s.',k,tStart,tStop);
    end
    segments(k).time = tk;
    segments(k).R = rk;
    segments(k).V = vk;
end

%% Compare truth and propagated ACD histories
truthACD = convert(t,r);
[trajectoryFig,ax] = makeFigure('Trajectory comparison','Position (m)');
[errorFig,ex] = makeFigure('Position error','Error (m)');
truthSeries = [truthACD,vecnorm(truthACD,2,2)];
for j = 1:4
    plot(ax(j),t,truthSeries(:,j),'k-','DisplayName','Truth');
end

edges = unique(windows(:));
for k = 1:numel(segments)
    ts = segments(k).time;
    rs = segments(k).R;
    propagatedACD = convert(ts,rs);
    % Match times in ECI before converting both histories at the same epoch.
    alignedTruthACD = convert(ts,interp1(t,r,ts,'linear'));
    delta = propagatedACD-alignedTruthACD;
    y = [propagatedACD,vecnorm(propagatedACD,2,2)];
    e = [delta,vecnorm(delta,2,2)];
    for j = 1:4
        p = plot(ax(j),ts,y(:,j),'--','Color',[0 .38 .72], ...
            'LineWidth',1.4,'DisplayName','Propagated');
        if k>1, p.HandleVisibility = 'off'; end
        plot(ex(j),ts,e(:,j),'Color',[0 .38 .72],'LineWidth',1.4);
        % Highlight the actual first error sample for each segment.
        plot(ex(j),ts(1),e(1,j),'o','MarkerSize',4,'Color',[0 .38 .72]);
    end
end
legend(ax(1),'Location','best');
for a = [ax ex]
    for k = 1:numel(edges)
        xline(a,edges(k),'k--','HandleVisibility','off');
    end
    xlim(a,edges([1 end]));
end
linkaxes(ax,'x');
linkaxes(ex,'x');

function [f,ax] = makeFigure(name,units)
f = figure('Name',name,'Color','w');
tl = tiledlayout(f,3,2,'TileSpacing','compact','Padding','compact');
ax = gobjects(1,4);
names = {'Altitude','Crossrange','Downrange','RSS of ACD components'};
for j = 1:3
    ax(j) = nexttile(tl,2*j-1);
end
ax(4) = nexttile(tl,2,[3 1]);
for j = 1:4
    hold(ax(j),'on');
    grid(ax(j),'on');
    title(ax(j),names{j});
    ylabel(ax(j),units);
end
xlabel(ax(3),'Time (s)');
xlabel(ax(4),'Time (s)');
title(tl,name);
end
