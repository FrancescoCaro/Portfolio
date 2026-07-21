# Homework 2 – Attitude Control System Design (BepiColombo @ Mercury)

Project developed for the **Space Missions and Systems** course (2026). The task is to design and verify a reaction-wheel-based Attitude Control System (ACS) that brings the BepiColombo spacecraft to a nadir-pointing attitude around Mercury and to assess how long the wheels can hold that attitude before saturating, under a solar radiation pressure (SRP) disturbance torque.

## Contents

| File | Description |
|---|---|
| `Homework2.pdf` | Original assignment text: orbital/geometric data, spacecraft and reaction-wheel parameters, and the three required tasks. |
| `Homework2_Caroletta.m` | Full MATLAB implementation solving all three tasks (single script, runs standalone). |
| `Caroletta_Homework_2_SMS.pdf` | Final report: derivations, methodology, results, plots and discussion. |

## Problem overview

The spacecraft orbits Mercury in a circular polar orbit and must reach nadir-pointing from a small initial attitude offset using 3 of its 4 tetrahedral reaction wheels, while respecting torque, current, power and settling-time constraints. Two disturbance scenarios are considered:

1. **Constant SRP torque** – solar panels always facing the Sun.
2. **Variable SRP torque** – solar panels fixed in the body frame, so the disturbance varies periodically with the orbital angle.

## What the script does

`Homework2_Caroletta.m` is organized in sections (MATLAB cell blocks, `%%`):

1. **Parameters** – orbital, planetary, SRP and reaction-wheel data from the assignment.
2. **Disturbance torque estimation** – computes the constant SRP torque `T_dy` from the geometry of the Sun/solar-panel/spacecraft vectors.
3. **Task 1 – Attitude control design** – integrates the closed-loop pitch dynamics (`ode113`) and tunes PD gains `Kp`, `Kd` via a `logspace()` search over candidate `Kp`, selecting the smallest gain pair that satisfies the ±4 arcsec / 5-minute pointing requirement and the wheel torque/current/power limits. Produces the 4-panel plot (pointing angle, wheel speeds, currents, powers).
4. **Task 2 – Desaturation time** – propagates wheel angular momentum analytically to estimate when Wheel 2/3 saturate, first ignoring eclipses and then accounting for the eclipse geometry (Mercury's shadow) to get a more realistic saturation time.
5. **Task 3 – Variable disturbance case** – redefines the SRP torque as a function of the attitude angle (`compute_Tdy`), shows that a PD controller alone cannot meet the requirements, and designs a **PID** controller (gains from `Kp/Kd/Ki` vs. natural frequency `ωn`) using the same iterative `logspace()` approach. Repeats the Task 1/2 analysis (pointing response and saturation time, including eclipses) for this case.

Helper functions defined inside the script: `dynamics`, `dynamics_new`, `dynamics_pid` (closed-loop ODEs for the three scenarios), `compute_Tdy` (SRP torque model), `saturation_event` (event function to stop integration at wheel saturation).

## How to run

Requires MATLAB (uses `ode113`, Symbolic-free). From the folder:

```matlab
Homework2_Caroletta
```

Running the script reproduces all numerical results (printed to the console) and generates all figures used in the report (pointing error, wheel telemetry, SRP torque profile, saturation analysis).

## Results summary

- A PD controller (`Kp ≈ 22.2`, `Kd ≈ 760`) meets the ±4 arcsec / 5-minute requirement under constant SRP disturbance; the wheels then saturate after **≈ 14.8 h** (accounting for eclipses).
- Under the variable (fixed-panel) disturbance, PD control alone is infeasible; a PID controller (`Kp ≈ 27.6`, `Kd ≈ 734`, `Ki ≈ 0.35`) is required, extending the saturation time to **≈ 18.6 h**.

Full derivations and discussion are in `Caroletta_Homework_2_SMS.pdf`.
