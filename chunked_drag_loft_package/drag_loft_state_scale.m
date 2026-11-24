function s = drag_loft_state_scale(Mach,Re,loftDeg)
% Install table and launch loft: drag_loft_state_scale('set',lut,loft)
% Query with current Mach and Re: drag_loft_state_scale(Mach,Re)
% Reset to nominal drag: drag_loft_state_scale('reset')
persistent model loft
if ischar(Mach) || isstring(Mach)
    switch lower(char(Mach))
        case 'set'
            if ~isscalar(loftDeg) || ~isfinite(loftDeg) || ...
                    loftDeg<min(Re.lofts) || loftDeg>max(Re.lofts)
                error('Loft outside calibrated range.');
            end
            model = Re;
            loft = loftDeg;
        case 'reset'
            model = []; loft = [];
        otherwise
            error('Unknown lookup command.');
    end
    s = []; return
end
if isempty(model), s = 1; return; end
p = ([loft,Mach,Re]-model.offset)./model.span;
s = model.F(p(1),p(2),p(3));
if ~isfinite(s) || s<=0
    error('DragLoft:OutsideCalibration', ...
        'No valid scale for loft %.6g deg, Mach %.6g, Re %.6g.',loft,Mach,Re);
end
end
