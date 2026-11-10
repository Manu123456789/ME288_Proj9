function [modeLookup, cases] = load_beta_inputs(cfg)
%LOAD_BETA_INPUTS Map production sources to the scaffold's input contract.
%
% TODO: load the validated tables from cfg.aeroSource. Prepare one lookup
% function per mode and return them in the cell vector modeLookup. Each
% function takes one query struct and returns the beta used by the existing
% drag equation. To replay the current single Mach table, capture its arrays:
% modeLookup = {@(q) interp1(Mtable, coefficientTable, q.Mach)};
%
% TODO: for the high-fidelity tables, prepare the actual axes and approved
% interpolation and domain rules. Each mode's function can use the fields
% supplied by the existing drag helper. Confirm the table value's units and
% normalization against the production force expression.
%
% TODO: map cfg.referenceSource to the cases struct array in the README.
% Include each case's initial conditions, times, and available truth signals.
% Use the scaffold units for comparison and apply any required frame mapping.

error('BetaScaffold:TODO', ...
    'Implement the table and truth mappings in load_beta_inputs.m. See README.md.');
end
