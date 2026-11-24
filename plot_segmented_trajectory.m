function [trajectoryFig,errorFig] = plot_segmented_trajectory(truth,segments,toACD)
% Compare recorded ECI histories over adjacent time windows.
% Example: [f1,f2] = plot_segmented_trajectory(truth,segments,@myACD);
% truth and each segments(k) have .time (absolute seconds on one clock)
% and .R (N-by-3 ECI position). Each segment includes both endpoints.
% Pass a function with this interface:
%   function acd = myACD(t,rECI)
%       % 1. Convert rECI at t using your ECI-to-ECEF routine.
%       % 2. Convert to [altitude,crossrange,downrange] in meters.
%       % Use the same epoch, origin, and range axes for every call.
%   end

[t,r] = history(truth);
assert(~isempty(segments),'Supply at least one recorded segment.');
truthACD = convert(toACD,t,r);
[trajectoryFig,ax] = makeFigure('Trajectory comparison','Position (m)');
[errorFig,ex] = makeFigure('Position error','Error (m)');
truthSeries = [truthACD,vecnorm(truthACD,2,2)];
for j = 1:4
    plot(ax(j),t,truthSeries(:,j),'k-','DisplayName','Truth');
end

edges = zeros(numel(segments)+1,1);
for k = 1:numel(segments)
    [ts,rs] = history(segments(k));
    assert(ts(1)>=t(1) && ts(end)<=t(end),'Segment exceeds truth coverage.');
    if k>1
        assert(abs(ts(1)-edges(k)) <= 16*eps(max(1,abs(ts(1)))), ...
            'Adjacent segments must share a boundary time.');
    end
    edges(k:k+1) = [ts(1);ts(end)];
    pred = convert(toACD,ts,rs);
    % Match times in ECI before converting both histories at the same epoch.
    ref = convert(toACD,ts,interp1(t,r,ts,'linear'));
    delta = pred-ref;
    y = [pred,vecnorm(pred,2,2)];
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
legend(ax(1),'show','Location','best');
for a = [ax ex]
    for k = 1:numel(edges)
        xline(a,edges(k),'k--','HandleVisibility','off');
    end
    xlim(a,edges([1 end]));
end
linkaxes(ax,'x');
linkaxes(ex,'x');
end

function [t,r] = history(s)
t = s.time(:);
r = s.R;
assert(isnumeric(t) && numel(t)>=2 && all(isfinite(t)) && all(diff(t)>0), ...
    'Each time vector must be finite and strictly increasing.');
assert(isnumeric(r) && isequal(size(r),[numel(t),3]) && all(isfinite(r(:))), ...
    'Each R must be a finite N-by-3 array.');
end

function acd = convert(toACD,t,r)
acd = toACD(t,r);
assert(isnumeric(acd) && isequal(size(acd),size(r)) && all(isfinite(acd(:))), ...
    'Conversion must return finite N-by-3 [altitude,crossrange,downrange].');
end

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
