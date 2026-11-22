function cL = lift_study_coefficient(t, knots, values)
%LIFT_STUDY_COEFFICIENT Fitted effective CL (independent of CFD tables).
% lift_study_coefficient('set', knots, values) installs a history.
% lift_study_coefficient(t) evaluates CL at time t.
% lift_study_coefficient('reset') disables fitted lift (CL = 0).

persistent timeNodes coeffNodes

if ischar(t) || isstring(t)
    if strcmpi(t, 'reset')
        timeNodes = []; coeffNodes = [];
    elseif strcmpi(t, 'set')
        [timeNodes, idx] = sort(knots(:));
        coeffNodes = values(idx);
    else
        error('Unknown lift coefficient command.');
    end
    cL = [];
    return
end

if isempty(coeffNodes)
    cL = 0;
else
    t = min(max(t, timeNodes(1)), timeNodes(end));
    cL = interp1(timeNodes, coeffNodes, t, 'linear');
end
end
