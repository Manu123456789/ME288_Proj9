function [a_drag, aero, cdMap] = Drag(R, V, t)
%DRAG Compute wind-axis drag acceleration for an ECI position and velocity.
%
%   a_drag = Drag(R, V)
%   [a_drag, aero] = Drag(R, V)
%   [a_drag, aero, cdMap] = Drag(R, V, t)
%
% R and V are three-element ECI vectors in SI units. Set cfg.alphaDeg to
% a scalar in degrees or a function of stage time t in seconds. The time
% input is required for a function. The first call caches the lookups.
% After changing settings or workspace data, run clear Drag.
% cdMap describes the Mach and Reynolds surface at the evaluated angle.
%
% Load CA_table and K_norm in the base workspace. Import the normal table:
%   CN_table = read_cn_table('your_CN_file.txt');
%   clear Drag
% CA_table is a struct vector with Reynolds, header, Mach, and data fields.
% CN_table can be one Reynolds-independent struct with header, Mach, data,
% or a legacy Reynolds-indexed struct vector. In both formats, column 1 of
% data is angle in degrees and the remaining columns match the Mach vector.
% Configure the Reynolds calculation below for the units of your axial table.
% At zero airspeed, a_drag is zero and aero is empty.

persistent model
if isempty(model)
    % User settings. Use '' for an optional table to use a zero coefficient.
    axialVar  = 'CA_table';
    normalVar = 'CN_table';
    sideVar   = '';             % Set to 'CY_table' to enable side force.

    cfg.alphaDeg = 0;                        % Constant angle of attack, deg
    % Fourier option. Use same-length rows for the cosine and sine terms.
    % tRef = 0;                           % s, in the propagator's clock
    % period = 8;                         % s
    % alphaMean = 4;                      % deg
    % aCos = [0.5, 0.2];                  % deg, one per harmonic
    % bSin = [1.5, 0];                    % deg, one per harmonic
    % omega = (2*pi/period) * (1:numel(aCos));
    % cfg.alphaDeg = @(t) alphaMean + sum(aCos .* cos(omega*(t - tRef)) ...
    %                                   + bSin .* sin(omega*(t - tRef)));

    cfg.K_norm = readDragVariable('K_norm');  % Existing S_ref/m in m^2/kg
    cfg.omegaEarth = [0, 0, 7e-5];            % Earth rate in rad/s
    cfg.altitudeFcn = @Altitude;             % H = f(R)
    cfg.atmosFcn = @Atmos;                   % [rho, a_sound] = f(H)
    cfg.reynoldsFcn = @reynoldsFromAtmos;     % Re = f(H, rho, U)
    cfg.bodyToWindFcn = @bodyToWindZeroSlip;  % [CD,CL,CS] = f(CA,CN,CY,alphaDeg,betaDeg)

    % To check a state with known Reynolds, temporarily set cfg.reynoldsFcn
    % to a function of (H, rho, U) that returns that in-range value.

    if ~isa(cfg.alphaDeg, 'function_handle')
        validateattributes(cfg.alphaDeg, {'numeric'}, {'scalar','real','finite'});
    end
    validateattributes(cfg.K_norm, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(cfg.omegaEarth, {'numeric'}, {'real','finite','numel',3});
    cfg.omegaEarth = reshape(cfg.omegaEarth, 1, 3);

    % Cache the lookup for each configured table.
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

alphaDeg = model.alphaDeg;
if isa(alphaDeg, 'function_handle')
    if nargin < 3
        error('Drag:TimeRequired', ...
            'Call Drag(R, V, t) when cfg.alphaDeg is a function of time.');
    end
    validateattributes(t, {'numeric'}, {'scalar','real','finite'});
    alphaDeg = model.alphaDeg(t);
    validateattributes(alphaDeg, {'numeric'}, {'scalar','real','finite'});
end

if nargout > 2
    cdMap.alphaDeg = alphaDeg;
    cdMap.tables = {model.CA, model.CN, model.CY};
    if ~isempty(model.CN) && isempty(model.CN.Re)
        % The existing plotter needs finite Re coverage for every surface.
        % Only the plot metadata repeats the same normal curve at both ends.
        reSpan = unique(model.CA.Re([1 end]));
        cdMap.tables{2}.Re = reSpan(:);
        cdMap.tables{2}.curves = repmat(model.CN.curves, numel(reSpan), 1);
    end
    cdMap.evaluate = @(Mach, Re) dragWindCoefficients(model, Mach, Re, alphaDeg);
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

[CD, CL, CS, CA, CN, CY] = dragWindCoefficients(model, Mach, Re, alphaDeg);
validateattributes(CD, {'numeric'}, {'scalar','real','finite','nonnegative'});

qbar = 0.5 * rho * U^2;
inverseBeta = model.K_norm * CD;
a_drag = reshape(-qbar * inverseBeta * (V_rel / U), size(V));

if nargout > 1
    aero = struct('Mach', Mach, 'Re', Re, 'alphaDeg', alphaDeg, ...
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
% Cache a Mach curve at a fixed angle or a Mach-by-angle grid per slice.

if ~isstruct(table) || isempty(table) || ~isvector(table)
    error('Drag:TableLayout', 'The %s table must be a nonempty struct vector.', coefficientName);
end
fields = fieldnames(table);
usesReynolds = numel(fields) >= 4;
if usesReynolds
    % Preserve the original positional field layout for Reynolds slices.
    reynoldsField = fields{1};
    machField = fields{3};
    dataField = fields{4};
    reynolds = zeros(numel(table), 1);
elseif strcmp(coefficientName, 'normal') && isscalar(table) && ...
        all(isfield(table, {'header','Mach','data'}))
    machField = 'Mach';
    dataField = 'data';
    reynolds = [];              % One angle/Mach grid, with no Reynolds axis.
else
    error('Drag:TableLayout', ...
        'Expected Reynolds/header/Mach/data slices or the scalar normal table from read_cn_table.');
end
machCurves = cell(numel(table), 1);
varyingAngle = isa(alphaDeg, 'function_handle');
for k = 1:numel(table)
    if usesReynolds
        reValue = table(k).(reynoldsField);
        validateattributes(reValue, {'numeric'}, {'scalar','real','finite','nonnegative'});
        reynolds(k) = reValue;
    end
    mach = table(k).(machField);
    data = table(k).(dataField);
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
    if ~varyingAngle && (alphaDeg < angles(1) || alphaDeg > angles(end))
        error('Drag:AngleOutsideTable', ...
            'Angle %g deg is outside the %s range [%g, %g] deg in slice %d.', ...
            alphaDeg, coefficientName, angles(1), angles(end), k);
    end

    coefficients = double(data(rowOrder, 1 + columnOrder));
    if varyingAngle
        % Keep Mach first so the CD map can still read GridVectors{1}.
        machCurves{k} = griddedInterpolant({mach, angles}, coefficients.', ...
            'linear', 'none');
    else
        atAngle = interp1(angles, coefficients, alphaDeg, 'linear');
        machCurves{k} = griddedInterpolant(mach, atAngle(:), 'linear', 'none');
    end
end

if usesReynolds
    [reynolds, order] = sort(reynolds);
    if any(diff(reynolds) <= 0)
        error('Drag:DuplicateReynolds', ...
            'The %s table must have unique Reynolds breakpoints.', coefficientName);
    end
    machCurves = machCurves(order);
end
lookup.Re = reynolds;
lookup.curves = machCurves;
lookup.coefficientName = coefficientName;
lookup.varyingAngle = varyingAngle;
end

function coefficient = dragCoefficientAt(lookup, Mach, Re, alphaDeg)
% A 2-D normal table uses its only curve. Reynolds-indexed tables blend slices.
% A fixed-angle lookup also accepts dragCoefficientAt(lookup, Mach, Re).

if isempty(lookup)
    coefficient = 0;
    return
end

reGrid = lookup.Re;
if isempty(reGrid)
    lower = 1;
    upper = 1;
else
    if Re < reGrid(1) || Re > reGrid(end)
        error('Drag:ReynoldsOutsideTable', ...
            'Re = %g is outside the %s range [%g, %g].', ...
            Re, lookup.coefficientName, reGrid(1), reGrid(end));
    end
    lower = find(reGrid <= Re, 1, 'last');
    upper = lower + (Re > reGrid(lower));
end
query = {Mach};
if lookup.varyingAngle
    query{2} = alphaDeg;
    for j = lower:upper
        curve = lookup.curves{j};
        angles = curve.GridVectors{2};
        if alphaDeg < angles(1) || alphaDeg > angles(end)
            error('Drag:AngleOutsideTable', ...
                'Angle %g deg is outside the %s range [%g, %g] deg in slice %d.', ...
                alphaDeg, lookup.coefficientName, angles(1), angles(end), j);
        end
    end
end

lowerCurve = lookup.curves{lower};
coefficient = lowerCurve(query{:});
if upper > lower
    upperCurve = lookup.curves{upper};
    fraction = (Re - reGrid(lower)) / (reGrid(upper) - reGrid(lower));
    coefficient = (1 - fraction) * coefficient + fraction * upperCurve(query{:});
end
if ~isfinite(coefficient)
    error('Drag:MachOutsideTable', ...
        'Mach %g is outside the required %s Mach grid(s).', Mach, lookup.coefficientName);
end
end

function [CD, CL, CS, CA, CN, CY] = dragWindCoefficients(model, Mach, Re, alphaDeg)
CA = dragCoefficientAt(model.CA, Mach, Re, alphaDeg);
CN = dragCoefficientAt(model.CN, Mach, Re, alphaDeg);
CY = dragCoefficientAt(model.CY, Mach, Re, alphaDeg);
[CD, CL, CS] = model.bodyToWindFcn(CA, CN, CY, alphaDeg, 0);
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
