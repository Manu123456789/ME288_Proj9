# 3-DOF trajectory and beta estimation scaffold

This MATLAB package connects validated empirical drag tables to an existing ECI trajectory propagator and compares its results against a collection of truth cases. The 3-DOF state consists of ECI position and velocity. The production propagator advances position with velocity and advances velocity with its existing gravity and axial drag accelerations. Table evaluation, timed mode selection, and validation are separate from those calculations.

## Files and integration work

| File | Status | Responsibility and required input |
|---|---|---|
| `beta_config.m` | Configuration template | Set table and reference sources, nominal step magnitude, shared switch times, and validation channels. Supply the beta unit and any agreed limits. |
| `load_beta_inputs.m` | Production mapping placeholder | Read your table and reference data, prepare each mode's lookup, and map truth histories and initial conditions to the contracts below. |
| `estimate_beta.m` | Complete scaffold | Select the active mode and call its prepared lookup to return scalar beta. |
| `propagate_beta_case.m` | Production connection placeholder | Connect to your existing propagator, pass the beta callback and schedule, and map the outputs to the history contract. |
| `run_beta_validation.m` | Complete scaffold | Call the propagator for each case, preserve case-specific inputs and statuses, and summarize errors across cases. |
| `compare_beta_case.m` | Complete scaffold | Align simulation outputs to truth times, compute residuals and metrics, and check configured peak limits. |
| `plot_beta_validation.m` | Complete scaffold | Plot ECI trajectory overlays, per-case error histories, and an empirical distribution of case RMS errors. |

The production handoff centers on `estimate_beta.m` and the table preparation and propagation connections you implement. The runner, comparisons, and plots support development and validation. The two placeholders raise a `BetaScaffold:TODO` error that identifies the missing connection.

## Run a validation batch

Place the files on the MATLAB path. Complete the table and truth mappings and connect the propagation adapter. Set `cfg.dt` and the shared or per-case switch times before running the batch.

```matlab
cfg = beta_config();
[modeLookup, cases] = load_beta_inputs(cfg);
results = run_beta_validation(cfg, modeLookup, cases);

metrics = struct2table(results.metrics);
spread = struct2table(rmfield(results.spread, 'values'));
plot_beta_validation(results, 'R');
```

For data already in memory, construct `modeLookup` and `cases` directly and pass them to the runner. You can also assign an existing adapter to `cfg.propagateCase` if it accepts the same arguments and returns the same history structure as `propagate_beta_case.m`.

## Drag model interface

Prepare `modeLookup` as a cell vector of function handles, one per flight mode, in chronological mode order. Each handle accepts a scalar `query` struct and returns one real, finite beta value in the table's confirmed convention. Choose the query field names, dimensions, units, interpolation method, and boundary behavior when connecting the actual tables. An evaluator can capture a prepared interpolant or call your existing lookup routine. Prepare data and interpolation structures once, before the force evaluation loop.

| Runtime input or output | Contract |
|---|---|
| `t` | Scalar time on the propagator's clock, in seconds. |
| `query` | Current table coordinates from the propagator's existing atmosphere and air-relative flow calculations. Include the actual axes, such as Mach, Reynolds number, angle, or an atmospheric quantity. |
| `modeLookup` | Prepared empirical lookup for each mode. The loader owns the mapping from the production data layout. |
| `switchTimes` | One strictly increasing finite event time for each transition between adjacent modes. |
| `segmentMode` | Optional table index held fixed during one integration segment. An empty value selects the mode by `t`. |
| `beta` | Scalar in the confirmed table normalization. |
| `mode` | Optional selected table index for diagnostics. |

In the development runner, `drag.betaAt(t, query, segmentMode)` binds the prepared lookup and the case schedule. Pass `[]` as the third argument for a point query. Direct calls to `estimate_beta` can omit the final `segmentMode` argument. Its schedule input must already follow the prepared contract.

For a confirmed rectilinear grid, MATLAB's `griddedInterpolant` is one way to prepare an evaluator over the actual axis vectors. The production lookup's interpolation and domain rules should guide that choice. A second mode can use a different query layout inside its own evaluator if the propagator supplies the needed fields. [MathWorks gridded interpolation reference](https://www.mathworks.com/help/matlab/ref/griddedinterpolant.html)

## Connect your existing propagator

`propagate_beta_case.m` receives the five inputs in your current interface and a sixth argument, `drag`. The first five are `R_0`, `V_0`, `t_0`, `t_final`, and `dt`. This adapter is where you translate the scaffold's column vectors, time conventions, units, and positive step magnitude to your production interface. The start and final times determine the propagation direction.

1. Make the new drag model available to the existing force evaluation. The proposed interface passes `drag` as an extra argument. Use the production integration point that fits your code.
2. Build `query` from the production ECI-to-flow calculations at each force evaluation, including intermediate integration stages. Use the existing atmosphere, Earth rotation, air-relative velocity, Mach, and any other quantities needed by the tables.
3. Replace the hard-coded coefficient lookup with `drag.betaAt(t, query, segmentMode)`. Use its scalar output in the confirmed existing drag equation, with the normalization already carried by the tables.
4. Use `drag.switchTimes` to place event boundaries in the propagation. Return the resulting state and acceleration histories in the common structure below. Log beta and flow channels where the production function supplies them.

### Table switching

With two tables, set one switch time. The point query convention selects Mode 1 before the event and Mode 2 at and after it. Additional modes use more tables and more ordered event times. The case can supply `cases(k).switchTimes` in place of `cfg.switchTimes`.

Split steps at switch times in the existing integrator. Hold the segment's mode for all stages, including the endpoint evaluation at the switch. A forward segment that ends at the first event therefore uses Mode 1 for that stage, and the next segment starts with Mode 2. Backward integration uses the same schedule and crosses from Mode 2 into Mode 1. Carry the propagated ECI state through the event and report the new mode for point diagnostics at the switch timestamp. Convert elapsed event times to the case's propagation clock during setup.

## Reference cases and histories

The case array uses an internal scaffold contract. Map production file names, variables, and data layouts in `load_beta_inputs.m`.

| Case field | Value |
|---|---|
| `id` | Unique nonempty character vector. |
| `R_0`, `V_0` | Three-element ECI vectors in m and m/s at `t_0`. |
| `t_0`, `t_final` | Requested start and end on the case clock, in s. |
| `truth` | Truth history in the structure below. |
| `switchTimes` | Optional case schedule. An empty field uses the shared schedule. |

Both `truth` and the adapter's returned `sim` contain `time` and `signals`. `time` is a finite vector of unique timestamps in s. Each field in `signals` has one row per timestamp, aligned with the supplied time order. Either time direction is supported. The comparator sorts samples internally.

| Signal | Shape | Scaffold convention |
|---|---|---|
| `R` | N × 3 | ECI position in m |
| `V` | N × 3 | ECI velocity in m/s |
| `acc_grav`, `acc_drag` | N × 3 each | ECI acceleration vectors in m/s² |
| `beta` | N × 1 | Confirmed scalar definition and unit, set in `cfg.channels` |
| `Re`, `Mach` | N × 1 each | Dimensionless values with agreed calculation conventions |

Include the channels available for each case. Omit absent fields and use `NaN` for missing samples within a present field. To add a quantity, add a row to `cfg.channels` with its field name, width, unit, required status, peak limit, and whether alignment should split at switch times. Supply comparable values in both histories. Convert both sources to the same unit and coordinate convention before comparing them.

## Validation and error spread

The runner calls the same prepared model for every case using that case's own initial conditions, end time, and schedule. It preserves case identifiers and reports propagation or comparison errors with their messages. Comparison uses truth samples within the requested window and linearly samples the model history there. Mode-dependent channels use samples from the same side of each switch. The alignment respects the available simulation time support and preserves exact matching samples. [MathWorks `interp1` reference](https://www.mathworks.com/help/matlab/ref/double.interp1.html)

`results.cases` holds both histories, signed component errors, error magnitudes, and the case status. `results.metrics` has one row per compared case and configured channel. Magnitudes are absolute error for scalar signals and Euclidean error for vectors. The metrics report sample RMS, peak, and the error at `t_final` when a matched truth sample exists there. `nMatched` and `nSamples` describe sample coverage.

A configured `peakLimit` is compared with the observed peak in the channel's unit. The margin is limit minus peak. A pass requires complete window and sample coverage. An observed exceedance reports a failure even with partial coverage. Required channels or channels with configured limits contribute to the case status. Unset limits leave available comparisons unassessed. Missing channels are reported individually, and beta metrics become available after its unit is set.

`results.spread` reports the mean, median, minimum, maximum, and individual values for each channel's RMS, peak, and endpoint metrics. Each complete case contributes once to a metric, and `nValid` and `nCases` show the population. The plot uses an empirical cumulative fraction for the selected channel's case RMS values. Error histories retain each case's elapsed propagation time. A pointwise mean or range across different time grids requires an agreed common comparison interval and alignment rule.

## Choice of beta model

| Form | When to use it |
|---|---|
| Direct lookup | Start here. The validated tables give beta as a function of current flight conditions and mode. Compare table nodes and representative interior queries against your reference implementation. |
| Fitted algebraic representation | Consider it when a simpler or faster representation is needed and measured lookup error and runtime justify the change. Use the tables as the fidelity reference. |
| Dynamic state-space representation | Consider it if data show repeatable dependence on prior conditions after accounting for current query inputs, mode, and known input errors. Establish the time basis and evidence for the proposed states first. |

The baseline beta map is algebraic. Its value changes as the existing propagator changes position, velocity, and the resulting flow conditions. A scheduled mode change selects a different map. A dynamic model would add a state vector, such as states representing an observed transient response, and would need a validated evolution law. In a local state-space formulation, `A` would govern those states, `B` would map deviations in agreed inputs into their evolution, `C` would map the states to beta deviations, and `D` would provide direct feedthrough. The operating point, initialization, mode transitions, and any state reset rules would need to be defined.

If a fitted representation is later justified, define a training set and independent validation cases. Scale coefficient and trajectory residuals by their respective units or tolerances before weighting them. Agree allowable deviations from the validated tables, assess the fitted result across operating conditions and both modes, and measure runtime in the production call path. The coefficient replay and trajectory comparisons give separate evidence for lookup fidelity and integrated accuracy.

## Outstanding engineering questions

### Critical before production integration

| Question | Information needed |
|---|---|
| What is beta? | Mathematical definition, units, sign, and the normalization already included in each table. |
| Which drag equation uses it? | The exact existing expression for drag acceleration, including dynamic pressure, beta's role, and the velocity direction convention. |
| How are the tables evaluated? | Actual axes and order per mode, units and reference quantities, approved interpolation method, and domain behavior. |
| Which frame and atmosphere conventions apply? | ECI definition, epoch, altitude, Earth rotation, wind conventions, air-relative velocity, density, sound speed, and Reynolds calculation. |
| Where does the model connect? | The actual propagation signature, force evaluation hook, time and step conventions, and available logged channels. |
| How are events scheduled? | The switch time source for each case, the time origin, exact boundary diagnostics, and the production solver's mechanism for splitting a crossing step. |
| Which truth data are available? | Each case's initial conditions, timestamps, frame and units, state and acceleration histories, and any beta, Mach, or Reynolds channels. Agree mandatory channels and evaluation windows before scoring compliance. |

### Decisions that can follow the direct lookup baseline

| Question | Evidence for the decision |
|---|---|
| Is a fitted algebraic representation needed? | Measured accuracy and runtime of the direct lookup in the existing propagator. |
| Is a state-space model justified? | Repeatable transient effects in beta under known inputs. Define states, inputs, outputs, `A`, `B`, `C`, `D`, operating points, and transition rules if that evidence exists. |
| How should fitting and performance assessment be organized? | Approved data splits, residual scales and weights, table-fidelity constraints, case acceptance limits, and a production runtime budget. |
| Is an across-case time band needed? | Actual sampling grids and durations, plus the agreed interval and alignment rule. |
