function [modeLookup, cases] = load_beta_inputs(cfg)
%LOAD_BETA_INPUTS Map production sources to the scaffold's input contract.
%
% TODO: load the validated tables from cfg.aeroSource. Prepare one lookup
% function per mode and return them in the cell vector modeLookup. Each
% function takes one query struct and returns a scalar beta. Confirm axes,
% interpolation, bounds, and normalization using the production definitions.
%
% TODO: map cfg.referenceSource to the cases struct array in the README.
% Include each case's initial conditions, times, and available truth signals.
% Use the scaffold units for comparison and apply any required frame mapping.

error('BetaScaffold:TODO', ...
    'Implement the table and truth mappings in load_beta_inputs.m. See README.md.');
end
