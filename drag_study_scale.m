function previous = drag_study_scale(value)
% Read the active drag multiplier, or set it and return the previous value.
persistent current
if isempty(current), current = 1; end
previous = current;
if nargin
    validateattributes(value, {'numeric'}, {'scalar','real','finite','positive'});
    current = value;
end
end
