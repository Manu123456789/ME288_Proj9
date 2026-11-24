# Chunked drag calibration and loft/Mach/Reynolds lookup

MATLAB R2021b syntax. The fitting study requires Optimization Toolbox (`lsqnonlin`).
No lift model is included. Your existing working Drag.m, Atmos, Altitude,
Reynolds calculation, aerodynamic tables, K_norm and five-input
propagate_traj (ode45) remain dependencies. Use one copy of each helper on
MATLAB's path.

## 1. Calibration (no change to Drag.m or propagate_traj)

Keep the original time-based multiplier active during the study:

```matlab
inverseBeta = drag_study_scale(t) * model.K_norm * CD;
```

Inputs: runs(j).label, loftDeg, time, R, V (Nx3 ECI position and velocity,
SI units, common physical time convention). The runs should already span the
intended calibration windows. For an original 30–150 s run:

```matlab
[out,chunkRuns] = fit_loft_drag_chunks(runs,dt,30,8);
```

The third input is chunk duration in seconds; the fourth is the initial
number of knots per chunk. Omit the fourth to use the original fitter's
default (20). The existing fitter still uses its own adaptive refinement
settings, regularization and 0.85–1.15 scale bounds, independently per chunk.
The splitter preserves the actual absolute timestamps, the original loft,
and all truth samples strictly inside each window. It interpolates BOTH R
and V at an exact boundary when that time is not present in the saved data.
The final chunk may be shorter. No propagated endpoint is carried forward.
This is a deterministic truth reset, not a Kalman filter measurement update.
The original fitter must be able to propagate to each requested endpoint
and return an output time vector covering the entire chunk.

`out` is a flat array with one element per chunk, not one per loft. Existing
fields (`label`, `loftDeg`, `time`, `scale`, `knots`, `knotScale`, `rmseR`,
`rmseV`) are preserved, and `parentCase`, `chunkIndex`, `interval` identify
the original run and window. `chunkRuns` has the matching truth histories.
Do not pass the original unsplit runs to the LUT builder after chunking.
Do not concatenate independent knot histories into a single time schedule.

## 2. Generate the state lookup (still in time-based mode)

```matlab
lut = build_drag_loft_lut(out,chunkRuns);
% Or associate fitted scales with corrected propagated states:
% lut = build_drag_loft_lut(out,chunkRuns,dt,'propagated');
save('drag_loft_chunked.mat','lut','out','chunkRuns');
```

The original builder and interpolator are unchanged. They read aero.Mach
and aero.Re from your Drag function, using its own air-relative speed and
Reynolds convention. The optional replay starts from truth for each chunk,
with the same model and solver settings as fitting. It repeats propagation,
not optimization. Truth sampling is the default; replay is useful when
corrected propagated states differ appreciably from truth.

Adjacent chunks share a truth endpoint, but their independent optimized
scales may disagree at that endpoint. The existing LUT builder averages
exact duplicate coordinates and reports `lut.duplicateSpread` (maximum
exact-duplicate disagreement); inspect this before treating the table as a
unique state-based function. `lut.samples.Case` indexes chunkRuns/out, so
original case and chunk identifiers remain recoverable. Repeated lofts
are permitted in the 3D scattered lookup. It requires a nondegenerate 3D
sample cloud spanning loft, Mach and Re; coverage is not assured simply
because the time interval is covered.

## 3. Runtime only

When switching from fitting to use of the already-generated lookup, change
only the existing multiplier line in Drag.m (same as the previous package):

```matlab
inverseBeta = drag_loft_state_scale(Mach,Re) * model.K_norm * CD;
```

Do not apply both corrections at once. Leave the rest of Drag.m and the
propagator unchanged. Switch back to the time-based line before further
fitting/replay. Clear Drag if needed to refresh your existing cached model.

```matlab
load('drag_loft_chunked.mat','lut');
drag_loft_state_scale('set',lut,16.7);
[R,V,time] = propagate_traj(R0,V0,t0,tf,dt);
drag_loft_state_scale('reset');
```

The state lookup retains the previous positive-scale and no-extrapolation
checks. A query outside the sample convex hull, including an internal
ode45 stage, throws an error. The hull itself does not establish model
accuracy. The persistent helpers are intended for serial simulations.

## 4. Validation

For a 30-second intended use, evaluate independently initialized windows
from held-out truth cases, including windows that straddle calibration
chunk boundaries. Compare to the uncorrected model and retain the original
continuous-duration calibration as a benchmark. Chunking reduces error
carryover in the fitting objective; it does not remove missing physics or
guarantee generalization. Continuous runs need separate validation if used.

Run `test_split_drag_runs` and `test_drag_loft_lut` for synthetic checks of
partitioning and interpolation. The files have not been executed in MATLAB
in this environment. Your actual fitting and runtime require your existing
propagator and aerodynamic data.
