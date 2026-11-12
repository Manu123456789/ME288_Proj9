function [a_drag, aero] = Drag(R, V)
%DRAG Compute wind-axis drag acceleration for an ECI position and velocity.
%
%   a_drag = Drag(R, V)
%   [a_drag, aero] = Drag(R, V)
%
% R and V are three-element ECI vectors in SI units. The first call reads
% the named base-workspace variables and prepares the fixed-angle lookup.
% The cached lookup is used for later calls. After changing a workspace
% table or K_norm, run clear Drag before the next propagation.
%
% Each table element contains Reynolds, header, Mach, and data fields in
% that order. Column 1 of data gives angle in degrees. The other columns
% give body-axis force coefficients in the order of the Mach vector.
% At zero airspeed, a_drag is zero and aero is empty.

persistent model
if isempty(model)
    % User settings. Use '' for an optional table to use a zero coefficient.
    axialVar  = 'CA_table';
    normalVar = 'CN_table';
    sideVar   = 'CY_table';

    cfg.alphaDeg = 0;
    cfg.K_norm = readDragVariable('K_norm');  % Existing S_ref/m in m^2/kg
    cfg.omegaEarth = [0, 0, 7e-5];            % Earth rate in rad/s
    cfg.altitudeFcn = @Altitude;             % H = f(R)
    cfg.atmosFcn = @Atmos;                   % [rho, a_sound] = f(H)
    cfg.reynoldsFcn = @reynoldsFromAtmos;     % Re = f(H, rho, U)
    cfg.bodyToWindFcn = @bodyToWindZeroSlip;  % [CD,CL,CS] = f(CA,CN,CY,alphaDeg,betaDeg)

    % To check a state with known Reynolds, temporarily set cfg.reynoldsFcn
    % to a function of (H, rho, U) that returns that in-range value.

    validateattributes(cfg.alphaDeg, {'numeric'}, {'scalar','real','finite'});
    validateattributes(cfg.K_norm, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(cfg.omegaEarth, {'numeric'}, {'real','finite','numel',3});
    cfg.omegaEarth = reshape(cfg.omegaEarth, 1, 3);

    % Prepare the angle slice once for each configured table.
    cfg.CA = prepareDragSlices(readDragVariable(axialVar), cfg.alphaDeg, 'axial');
    cfg.CN = [];
    if ~isempty(normalVar)
        cfg.CN = prepareDragSlices(readDragVariable(normalVar), cfg.alphaDeg, 'normal');
    end
    cfg.CY = [];
    if ~isempty(sideVar)
        cfg.CY = prepareDragSlices(readDragVariable(sideVar), cfg.alphaDeg, 'side');
    end
    model = cfg;
end

% Air-relative velocity in ECI.
V_atm = cross(model.omegaEarth, reshape(R, 1, 3));
V_rel = reshape(V, 1, 3) - V_atm;
U = norm(V_rel);
if U == 0
    a_drag = zeros(size(V));
    aero = [];
    return
end

H = model.altitudeFcn(R);
[rho, a_sound] = model.atmosFcn(H);
validateattributes(rho, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(a_sound, {'numeric'}, {'scalar','real','finite','positive'});

Mach = U / a_sound;
Re = model.reynoldsFcn(H, rho, U);
validateattributes(Re, {'numeric'}, {'scalar','real','finite','nonnegative'});

CA = dragCoefficientAt(model.CA, Mach, Re);
CN = dragCoefficientAt(model.CN, Mach, Re);
CY = dragCoefficientAt(model.CY, Mach, Re);
[CD, CL, CS] = model.bodyToWindFcn(CA, CN, CY, model.alphaDeg, 0);
validateattributes(CD, {'numeric'}, {'scalar','real','finite','nonnegative'});

qbar = 0.5 * rho * U^2;
inverseBeta = model.K_norm * CD;
a_drag = reshape(-qbar * inverseBeta * (V_rel / U), size(V));

if nargout > 1
    aero = struct('Mach', Mach, 'Re', Re, 'alphaDeg', model.alphaDeg, ...
        'CA', CA, 'CN', CN, 'CY', CY, 'CD', CD, 'CL', CL, 'CS', CS, ...
        'qbar', qbar, 'inverseBeta', inverseBeta);
end
end

function value = readDragVariable(name)
% Read a named variable from the MATLAB base workspace.
if ~isvarname(name) || evalin('base', sprintf('exist(''%s'', ''var'')', name)) ~= 1
    error('Drag:MissingWorkspaceVariable', ...
        'Load the variable %s into the base workspace before calling Drag.', name);
end
value = evalin('base', name);
end

function lookup = prepareDragSlices(table, alphaDeg, coefficientName)
% Interpolate the angle rows once. Keep each Reynolds slice's Mach grid.

if ~isstruct(table) || isempty(table) || ~isvector(table) || numel(fieldnames(table)) < 4
    error('Drag:TableLayout', ...
        'The %s table must be a struct vector with Reynolds, header, Mach, and data fields.', ...
        coefficientName);
end

fields = fieldnames(table);
reynoldsField = fields{1};
machField = fields{3};
dataField = fields{4};

reynolds = zeros(numel(table), 1);
machCurves = cell(numel(table), 1);
for k = 1:numel(table)
    reValue = table(k).(reynoldsField);
    mach = table(k).(machField);
    data = table(k).(dataField);
    validateattributes(reValue, {'numeric'}, {'scalar','real','finite','nonnegative'});
    validateattributes(mach, {'numeric'}, {'vector','real','finite','nonnegative'});
    validateattributes(data, {'numeric'}, {'2d','real','finite'});

    if numel(mach) < 2 || size(data, 1) < 2 || size(data, 2) ~= numel(mach) + 1
        error('Drag:SliceDimensions', ...
            '%s slice %d needs angle rows and one coefficient column per Mach value.', ...
            coefficientName, k);
    end

    [angles, rowOrder] = sort(double(data(:,1)));
    [mach, columnOrder] = sort(double(mach(:)));
    if any(diff(angles) <= 0) || any(diff(mach) <= 0)
        error('Drag:DuplicateBreakpoints', ...
            '%s slice %d must have unique angle and Mach breakpoints.', coefficientName, k);
    end
    if alphaDeg < angles(1) || alphaDeg > angles(end)
        error('Drag:AngleOutsideTable', ...
            'Angle %g deg is outside the %s slice at Re = %g.', ...
            alphaDeg, coefficientName, reValue);
    end

    coefficients = double(data(rowOrder, 1 + columnOrder));
    atAngle = interp1(angles, coefficients, alphaDeg, 'linear');
    machCurves{k} = griddedInterpolant(mach, atAngle(:), 'linear', 'none');
    reynolds(k) = reValue;
end

[lookup.Re, order] = sort(reynolds);
if any(diff(lookup.Re) <= 0)
    error('Drag:DuplicateReynolds', ...
        'The %s table must have unique Reynolds breakpoints.', coefficientName);
end
lookup.curves = machCurves(order);
lookup.coefficientName = coefficientName;
end

function coefficient = dragCoefficientAt(lookup, Mach, Re)
% Evaluate Mach on the bracketing slices, then interpolate in Reynolds.

if isempty(lookup)
    coefficient = 0;
    return
end

reGrid = lookup.Re;
if Re < reGrid(1) || Re > reGrid(end)
    error('Drag:ReynoldsOutsideTable', ...
        'Re = %g is outside the %s range [%g, %g].', ...
        Re, lookup.coefficientName, reGrid(1), reGrid(end));
end

lower = find(reGrid <= Re, 1, 'last');
lowerCurve = lookup.curves{lower};
coefficient = lowerCurve(Mach);
if Re > reGrid(lower)
    upperCurve = lookup.curves{lower + 1};
    fraction = (Re - reGrid(lower)) / (reGrid(lower + 1) - reGrid(lower));
    coefficient = (1 - fraction) * coefficient + fraction * upperCurve(Mach);
end
if ~isfinite(coefficient)
    error('Drag:MachOutsideTable', ...
        'Mach %g at Re = %g is outside the required %s Mach grid(s).', ...
        Mach, Re, lookup.coefficientName);
end
end

function [CD, CL, CS] = bodyToWindZeroSlip(CA, CN, CY, alphaDeg, ~)
% Positive axial aft, normal up, side right. Called at zero sideslip.
CD = CA * cosd(alphaDeg) + CN * sind(alphaDeg);
CL = CN * cosd(alphaDeg) - CA * sind(alphaDeg);
CS = CY;
end

function Re = reynoldsFromAtmos(H, rho, U) %#ok<INUSD,STOUT>
% REYNOLDS CALCULATION TO FILL IN
% H is altitude, rho is density in kg/m^3, and U is air-relative speed in m/s.
% With mu in Pa*s and the table reference length L_ref in m, use
%     Re = rho * U * L_ref / mu;
% For a unit Reynolds table, use rho * U / mu in the table's inverse-length units.
error('Drag:ReynoldsNotConfigured', ...
    'Fill in reynoldsFromAtmos or set cfg.reynoldsFcn in Drag.m.');
end
