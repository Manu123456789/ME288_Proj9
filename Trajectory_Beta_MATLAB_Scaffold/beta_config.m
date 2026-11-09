function cfg = beta_config()
%BETA_CONFIG Settings for data mapping, propagation, and validation.

cfg.aeroSource = [];       % TODO: locations or in-memory inputs for the two tables.
cfg.referenceSource = [];  % TODO: truth histories and their initial conditions.
cfg.dt = [];               % TODO: positive nominal step magnitude in seconds.
cfg.switchTimes = [];      % TODO: one event time per mode transition, in seconds.
cfg.propagateCase = @propagate_beta_case;

% Case-specific switchTimes, when supplied, take precedence over this schedule.
% TODO: set the beta unit after agreeing its meaning and normalization.
% TODO: set peakLimit in each channel's unit when a limit is agreed.
cfg.channels = struct( ...
    'name',          {'R', 'V',   'acc_grav', 'acc_drag', 'beta', 'Re', 'Mach'}, ...
    'width',         {3,   3,     3,          3,          1,      1,    1}, ...
    'unit',          {'m', 'm/s', 'm/s^2',    'm/s^2',    '',     '1',  '1'}, ...
    'required',      {true,true,  false,      false,      false,  false,false}, ...
    'peakLimit',     {[],  [],    [],         [],         [],     [],   []}, ...
    'splitAtSwitch', {false,false,false,     true,       true,   true, false});
end
