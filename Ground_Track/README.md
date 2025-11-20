# Satellite Ground Track Plotter (MATLAB)

A simple, lightweight, and fully commented MATLAB toolbox to compute and visualize the **ground track** of a satellite starting from classical orbital elements.

## Features
- Converts classical orbital elements → ECI position & velocity  
- High-precision two-body propagation using `ode45`  
- Accurate sub-satellite point calculation (accounts for Earth rotation)  
- Plots the ground track on a high-resolution Mercator world map  
- Marks the initial point in green  
- Optional diagnostic plots (latitude, ECI/ECEF longitude vs time)

## Requirements
- MATLAB R2016b or newer  
- Image Processing Toolbox (only for `imread`/`imagesc`)  
- No external dependencies

## Files
| File                       | Description                                                                 |
|----------------------------|-----------------------------------------------------------------------------|
| `COEtor0v0.m`              | Converts COE (a, e, i, Ω, ω, ν) → initial ECI position & velocity [km, km/s] |
| `dyn_mod.m`                | Two-body dynamics function for `ode45` (Keplerian only)                     |
| `tracciaaterra.m`          | Main script – propagation, ground-track calculation and plotting           |
| `Mercator_projection_SW.*` | Equirectangular world map (jpg/png/bmp – filename must start with `Mercator_projection_SW`) |

## Quick Start

1. Place all files in the same folder.  
2. Put your Mercator map in the same folder (e.g. `Mercator_projection_SW.jpg`).  
3. Open `tracciaaterra.m` and set your orbital parameters (top of the file):

```matlab
% Example: ISS-like orbit
a    = 6378 + 420;   % semi-major axis [km]
ecc  = 0.0005;       % eccentricity
i    = 51.6;         % inclination [deg]
raan = 120;          % RAAN [deg]
w    = 90;           % argument of perigee [deg]
v    = 0;            % true anomaly [deg] (0 = perigee)

[r0_COE, v0_COE] = COEtor0v0(a, ecc, i, raan, w, v);
```
4. Run `tracciaaterra.m`

You will get: 
  * Ground track on the world map (red dots, green start point)
  * Three diagnostic subplots
