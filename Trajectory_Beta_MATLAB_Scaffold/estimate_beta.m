function [beta, mode] = estimate_beta(t, query, modeLookup, switchTimes, segmentMode)
%ESTIMATE_BETA Evaluate a prepared empirical lookup for the active mode.
% switchTimes is the prepared, strictly increasing schedule on the t clock.
% query has the fields expected by the production lookup functions.
% The current interp1 baseline uses query.Mach. A high-fidelity table can
% consume additional coordinates through its prepared lookup function.
% segmentMode selects the table for all stages of an integration segment.
% Supply [] or omit segmentMode for a time-based point query.

validateattributes(t, {'numeric'}, {'real', 'finite', 'scalar'});
if nargin < 5 || isempty(segmentMode)
    mode = 1 + nnz(t >= switchTimes);
else
    mode = segmentMode;
end
validateattributes(mode, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', '>=', 1, '<=', numel(modeLookup)});

lookup = modeLookup{mode};
beta = lookup(query);
validateattributes(beta, {'numeric'}, {'real', 'finite', 'scalar'}, ...
    mfilename, 'beta returned by the mode lookup');
end
