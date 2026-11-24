function s = drag_study_scale(t, knots, values)
%DRAG_STUDY_SCALE Time-varying multiplier for Drag.m.
% drag_study_scale('set', timeNodes, scaleNodes) installs a schedule.
% drag_study_scale(t) returns the multiplier at time t.
% drag_study_scale('reset') restores scale = 1.

persistent timeNodes scaleNodes

if ischar(t) || isstring(t)
    if strcmpi(t, 'reset')
        timeNodes = []; scaleNodes = [];
    elseif strcmpi(t, 'set')
        [timeNodes, idx] = sort(knots(:));
        scaleNodes = values(idx);
    else
        error('Unknown drag scale command.');
    end
    s = [];
    return
end

if isempty(scaleNodes)
    s = 1;
else
    t = min(max(t, timeNodes(1)), timeNodes(end));
    s = interp1(timeNodes, scaleNodes, t, 'linear');
end
end
