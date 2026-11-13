# Drag coefficient map and trajectory queries

The map shows Mach on the horizontal axis, Reynolds on the vertical axis,
and wind-axis drag coefficient `CD` in color. Each red dot represents a
saved trajectory state. The connected line and arrowheads follow the
order of the supplied state rows. Click a red dot for its time, Mach,
Reynolds coordinate, and `CD`.

## Add it to your package

1. Put `plot_drag_cd_map.m` on your MATLAB path. The data tips use MATLAB
   R2019a or newer.
2. Add the optional third output to your existing `Drag.m` using the
   three edits below. The `reference/Drag.m` file and
   `Drag_map_changes.diff` show the edits against the original.
   Apply them to the copy that has your atmosphere and Reynolds calculation.
3. Run the five-input propagator, then pass its saved states to the plot.

```matlab
[R,V,time,acc_grav,acc_drag] = propagate_traj(R_o,V_o,t0,final_time,dt);
map = plot_drag_cd_map(time, R, V);
```

The plot expects one time value per row of the `N×3` ECI position and
velocity histories. Pass the rows in the order the propagator produced
them. For a backward run, the arrows follow decreasing time.

If the table stores raw Reynolds values, the following example displays
millions and shows every fifth red dot plus the last valid dot. The
path, arrow directions, and returned data use every saved state.

```matlab
map = plot_drag_cd_map(time, R, V, @Drag, 1e6, 5);
```

Use a divisor of 1 when the table values are already in the scale you
want to show. The divisor changes the display only.

### Optional third output in Drag

Change the first line.

```matlab
function [a_drag, aero, cdMap] = Drag(R, V)
```

Immediately after the initialization block and before the airspeed check,
insert the following. This makes the surface available even when the first
saved state has zero relative airspeed.

```matlab
if nargout > 2
    cdMap.alphaDeg = model.alphaDeg;
    cdMap.tables = {model.CA, model.CN, model.CY};
    cdMap.evaluate = @(Mach, Re) dragWindCoefficients(model, Mach, Re);
end
```

Replace the three coefficient lookups and the body-to-wind call in the
main function with the shared call below. Keep the existing `CD` validation
and acceleration calculation after it.

```matlab
[CD, CL, CS, CA, CN, CY] = dragWindCoefficients(model, Mach, Re);
```

Add this local function at the end of the file. If you use a custom
coefficient conversion, move your existing calculation into this helper.
The map needs a converter that depends on the coefficients and fixed angles
in this lookup.

```matlab
function [CD, CL, CS, CA, CN, CY] = dragWindCoefficients(model, Mach, Re)
CA = dragCoefficientAt(model.CA, Mach, Re);
CN = dragCoefficientAt(model.CN, Mach, Re);
CY = dragCoefficientAt(model.CY, Mach, Re);
[CD, CL, CS] = model.bodyToWindFcn(CA, CN, CY, model.alphaDeg, 0);
end
```

If `Drag` is local to `propagate_traj.m`, place the following line inside
`propagate_traj` after the histories have been computed. `@Drag` then
resolves to the same local function used for propagation.

```matlab
plot_drag_cd_map(time, R, V, @Drag);
```

## Read the plot

The surface samples the same fixed-angle lookup and wind-axis conversion
as drag. At zero angle of attack and zero sideslip, `CD = CA` with the
default converter. When adjacent Reynolds slices have different Mach
grids, a colored strip appears only where all needed curves cover that
query. Exact Reynolds breakpoints are also shown. Keep the same cached
table configuration for the propagation and replay so both use the same
coefficients. If you change the tables, rebuild the cache and rerun the
trajectory before plotting.

The line connects the saved-state queries in order, and the arrowheads
indicate direction along the line. The plot re-evaluates drag at the saved
position and velocity rows after the run. Intermediate Heun evaluations
require their own logged states. A zero-airspeed row gives a gap in the
path because the lookup is skipped at that state.

`map.time`, `map.Mach`, `map.Re`, and `map.CD` give the results aligned to
the input history. Click a dot for the exact values, or rotate the MATLAB
axes to see the height of the surface as well as its color.
For example, use `view(map.axes, 42, 28)`.
`Drag_CD_Map_Example.png` illustrates the display with synthetic data.
