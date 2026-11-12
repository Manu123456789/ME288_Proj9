# Drag.m README

`Drag(R,V)` takes ECI position in meters and velocity in m/s. It returns ECI drag acceleration in m/s² with the same row or column shape as `V`. Heun calls it at the current and predicted states. The earlier `makeDragFcn` returned a callable handle. This version prepares the lookup on the first `Drag` call and keeps the setup inside `Drag.m`.

## Run it

1. Put `Drag.m` on the MATLAB path. If `propagate_traj.m` defines a local `Drag`, replace that definition and include the local helpers from `Drag.m`.
2. Load `CA_table` and `K_norm` into the base workspace. Load `CN_table` and `CY_table` for the enabled normal and side tables, or set `normalVar = ''` and `sideVar = ''` for axial-only use. Set the variable names and fixed angle at the top of `Drag.m`.
3. Fill in **`reynoldsFromAtmos(H,rho,U)`** at the bottom. Then use your five-input propagator call.

```matlab
[R,V,time,acc_grav,acc_drag] = propagate_traj(R_o,V_o,t0,final_time,dt);
```

## Table layout

A `1×24` struct has 24 entries. Each entry represents a Reynolds slice and carries four fields. The reader also accepts column struct arrays and different numbers of slices in each coefficient table.

| Field position | Contents |
|---|---|
| 1 | The Reynolds value for that slice, in the same convention as the runtime query |
| 2 | Header text |
| 3 | Mach breakpoints, typically a `1×19` vector |
| 4 | Typically a `39×20` matrix. Column 1 holds angles in degrees. Columns 2–20 hold coefficients in the order of the Mach vector. |

Run `whos CA_table`, `fieldnames(CA_table)`, and `CA_table(1)` to check the layout. If the first stored field contains record IDs such as 1, 2, 3, identify the Reynolds value for each record before using the lookup. The fixed angle must be within every configured slice. Runtime Mach and Reynolds queries must be in range for each enabled table.

## Settings

| Setting | Purpose |
|---|---|
| `cfg.alphaDeg = 0` | Fixed angle of attack in degrees. The same angle selects the lookup slice and goes into the force rotation. |
| `cfg.K_norm` | Reference area divided by mass, in m²/kg. Use the area that defines the table coefficients and the mass used by the propagator. For the old `qbar/beta` form, `K_norm * CD` gives `1/beta`. |
| `cfg.omegaEarth` | Earth rotation vector in rad/s. Use your propagator's Earth rate for the co-rotating atmosphere model. |
| `cfg.altitudeFcn = @Altitude` | Calls your altitude function with ECI position `R` to get `H`. |
| `cfg.atmosFcn = @Atmos` | Calls your atmosphere function with `H`. Its first two outputs must be density `rho` and sound speed `a_sound`, in that order. |
| `cfg.reynoldsFcn` | Calls the Reynolds calculation with `H`, `rho`, and air-relative speed `U`. |
| `cfg.bodyToWindFcn` | Uses the included zero-sideslip conversion by default. A custom converter takes `(CA,CN,CY,alphaDeg,betaDeg)` and returns `[CD,CL,CS]`. |

The `@` creates a function handle so the settings can point to your existing routines. Any gas-property gamma used to compute sound speed is handled by the atmosphere routine, which supplies `a_sound`.

For a length-based Reynolds grid, use `Re = rho*U*L_ref/mu` with the table's reference length and dynamic viscosity. For a unit-Reynolds grid, use `rho*U/mu` in the table's inverse-length units. A known in-range Reynolds value can be used temporarily to check a single state.

## What Drag does

On the first call, the code reads the workspace inputs and prepares one Mach curve at `alphaDeg` for each Reynolds slice. Each drag evaluation subtracts `cross(omegaEarth,R)` from `V` to get air-relative velocity. It uses that vector's speed `U` to get `Mach = U/a_sound`, Reynolds number, and `qbar = 0.5*rho*U^2`. The lookup interpolates in Mach on the bracketing slices and blends the results in Reynolds. The wind-axis `CD` sets the acceleration magnitude `qbar * K_norm * CD` opposite the air-relative velocity.

The included rotation assumes positive axial aft, normal up, and side right. At zero sideslip it gives `CD = CA*cosd(alphaDeg) + CN*sind(alphaDeg)` and `CS = CY`. At zero angle, `CA` alone supplies drag. At nonzero angle, `CN` contributes to drag. `CL` and `CS` are diagnostics.

With the standalone file, use `[a_drag,aero] = Drag(R_i,V_i)` to inspect Mach, Reynolds, coefficients, dynamic pressure, and inverse ballistic coefficient for one state. At zero airspeed, acceleration is zero and `aero` is empty. After changing workspace tables or `K_norm`, run `clear Drag` before the next run. Use `clear propagate_traj` if `Drag` is local to that file.
