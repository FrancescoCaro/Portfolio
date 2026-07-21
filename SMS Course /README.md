# Space Missions and Systems (2026)

This repository collects the assignments completed for the **Space Missions and Systems** course. Each assignment is a self-contained project (MATLAB code + PDF report) exploring a different aspect of spacecraft navigation and control around Mercury, using the **BepiColombo / Mercury Planetary Orbiter (MPO)** mission as the reference scenario. Every subfolder has its own detailed README; this page is just an index.

## Assignments

| # | Folder | Topic | Summary |
|---|---|---|---|
| 1 | `Orbit Determination/` | Orbit Determination | Batch least-squares estimation of the MPO state (position, velocity, Mercury's GM, solar radiation pressure coefficient) around Mercury from simulated range and range-rate tracking data, including uncertainty propagation and analysis of a data gap due to occultation. |
| 2 | `Attitude Control/` | Attitude Determination and Control | Design of a reaction-wheel attitude control system to achieve nadir-pointing under solar radiation pressure disturbance, including PD/PID gain tuning and reaction wheel desaturation/saturation-time analysis. |
| 3 | `Challenge/` | Advanced Orbit Determination (Challenge) | Extended orbit determination exercise for an Earth orbiter: batch least-squares with atmospheric drag estimation, followed by an Extended Kalman Filter (EKF) using GNSS pseudo-range/pseudo-range-rate observables to also estimate the onboard clock bias and drift. |

## Course context

The course covers the mission analysis, navigation and control aspects of spacecraft operations, using BepiColombo's cruise and Mercury orbit-insertion phases as a running case study. Assignments progress from orbit determination fundamentals (Homework 1) to attitude control system design (Homework 2), and finally to a more advanced orbit determination challenge combining batch estimation with sequential (EKF) filtering (Challenge).

## Structure and conventions

Each assignment folder follows the same structure:
- **Assignment text (PDF)** – original problem statement.
- **MATLAB code (`.m`)** – self-contained script, runs with a single command in MATLAB, prints result tables to the console and generates all report figures.
- **Report (PDF)** – full write-up with derivations, methodology, results and discussion.

All simulations use `ode113` with high-precision tolerances (`RelTol = AbsTol = 1e-13`), as recommended for orbital dynamics problems.
