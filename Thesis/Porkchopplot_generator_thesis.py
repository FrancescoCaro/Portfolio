import numpy as np
import spiceypy as spice
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import pykep as pk
from pykep import propagate_lagrangian

# Importing the NASA/NAIF kernels 
spice.furnsh("/mnt/c/Users/Francesco/Desktop/Tesi/Gravity.tpc")  # Contains the values of GM
spice.furnsh("/mnt/c/Users/Francesco/Desktop/Tesi/naif0012.tls")  # LSK (Leap Seconds Kernel)
spice.furnsh("/mnt/c/Users/Francesco/Desktop/Tesi/de440s.bsp")    # SPK (Planet Ephemeris Kernel)

# Let us define C₃, v∞, and a Δv cut-off threshold beyond which we are no longer interested in performing the calculations
cutoff_c3 = 50
cutoff_dv = 20
cutoff_v_inf =  12

### DEFINITION OF THE LAUNCH AND ARRIVAL DATE RANGES
# Definition of the launch window
launch_start_day =  "2020-01-01 00:00:00 UTC" 
launch_end_day  = "2020-12-31 00:00:00 UTC"

# Definition of the arrival window
arrival_start_day = "2020-10-01 00:00:00 UTC"
arrival_end_day  = "2021-12-12 00:00:00 UTC"

step = 1 # Time step
days2sec = 86400.0
step_s = step*days2sec

# Conversion of the window boundaries to Ephemeris Time (ET) (seconds past J2000)
et_launch_start = spice.utc2et(launch_start_day)
et_launch_end   = spice.utc2et(launch_end_day)
et_arrival_start = spice.utc2et(arrival_start_day)
et_arrival_end   = spice.utc2et(arrival_end_day)

# Construction of the launch and arrival date vectors
launch_ets = np.arange(et_launch_start, et_launch_end +1, step_s, dtype=float)
arrival_ets = np.arange(et_arrival_start, et_arrival_end +1, step_s, dtype=float)

# Number of days for each window and total combinations
ls = len(launch_ets)
as_ = len(arrival_ets)
total_c = ls*as_

### ALGORITHM IMPLEMENTATION

# Load the gravitational constants [km^3]/[s^2]
mu_sun   = spice.gdpool("BODY10_GM", 0, 1)[0]
mu_earth = spice.gdpool("BODY399_GM", 0, 1)[0]
mu_mars  = spice.gdpool("BODY499_GM", 0, 1)[0] 

# Computation of the relevant orbital parameters for Δv calculation
r_LEO = (6371 + 300)      # Earth radius + LEO 300 km [km]
v_esc = np.sqrt(2*mu_earth/r_LEO)        # escape velocity from a circular orbit of radius r_LEO [km/s]

r_LMO = (3390 + 200) # radius of the final circular orbit of Mars [km]
v_esc_M =  np.sqrt(2*mu_mars/r_LMO)  # escape velocity from a Mars circular orbit of radius r_LMO [km/s]


# Initialize the vectors containing the velocities and TOFs
C3 = np.full((as_,ls),np.nan)    # Matrix with the values of C3
tofs = np.full((as_,ls),np.nan)   # Matrix with the values TOFs
Dv = np.full((as_,ls),np.nan)     # Matrix with the values of Dv
v_infs = np.full((as_,ls),np.nan)  # Matrix with the values of  v_inf_arr

## Main loop to build matrices with TOF, C3, v_inf_arr and DeltaV
for i, launch_date in enumerate(launch_ets):

   # Obtain the planet positions in the heliocentric frame
    state_earth, lt = spice.spkezr('EARTH BARYCENTER', launch_date, 'ECLIPJ2000', 'NONE', 'SUN') 
    r_earth = state_earth[:3]   # Position vector of the Earth in the Solar reference frame [km]
    v_earth = state_earth[3:]   # Velocity vector of the Earth in the solar reference frame [km/s]

    for j, arrival_date in enumerate(arrival_ets):

        # Skip combinations where the arrival date is earlier than the launch date
        if arrival_date <= launch_date:
            continue

        # Compute the values of TOF for each launch_date and arrival_date in seconds [s]
        TOF = arrival_date - launch_date  
        tofs[j, i] = TOF  # Store the value of  TOF in the matrix tofs

        state_mars, lt = spice.spkezr('MARS BARYCENTER', arrival_date, 'ECLIPJ2000', 'NONE', 'SUN')
        r_mars = state_mars[:3]     # Position vector of Mars in the Solar reference frame [km]
        v_mars = state_mars[3:]     # Velocity vector of Mars in the solar reference frame [km/s]

        try:
            # Solve the Lambert problem between Earth and Mars for the computed TOF
            sol = pk.lambert_problem(r_earth, r_mars, TOF, mu_sun)
            # Extract the departure and arrival velocities from the Lambert solution
            v1 = np.array(sol.get_v1()[0])  # departure velocity from Earth [km/s]
            v2 = np.array(sol.get_v2()[0])  # arrival velocity at Mars [km/s]
        except:
            continue

        # Compute C₃ (departure) and v∞ (arrival)
        c3 = (np.linalg.norm(v1 - v_earth))**2  # [km^2/s^2]
        v_inf_arr = np.linalg.norm(v2 - v_mars)  # [km/s]

        # Computation of the required Δv
        dv1 = np.sqrt(c3 + v_esc**2) - np.sqrt(mu_earth/r_LEO)
        dv2 = np.sqrt(v_inf_arr**2 + v_esc_M**2) - np.sqrt(mu_mars/r_LMO)
        dv = dv1 + dv2  # Total DeltaV [km/s]

        # Apply the C₃, v∞_arr and Δv cutoffs
        c3 = min(c3, cutoff_c3)
        v_inf_arr = min(v_inf_arr, cutoff_v_inf)
        dv = min(dv, cutoff_dv)

        # Store the values in the corresponding matrices
        C3[j, i] = c3   
        v_infs[j, i] = v_inf_arr  
        Dv[j, i] = dv

# Conversion ET → datetime 
launch_dates = [datetime.fromisoformat(spice.et2utc(et, 'ISOC', 0)) for et in launch_ets]
arrival_dates = [datetime.fromisoformat(spice.et2utc(et, 'ISOC', 0)) for et in arrival_ets]

### PLOT THE RESULTS

# Grid coordinates creation 
X, Y = np.meshgrid(launch_dates, arrival_dates)

#  Define contour levels for the plots 
c3_levels = np.arange(0, cutoff_c3, 5)
v_inf_levels = np.arange(0, cutoff_v_inf, 3)
tof_levels_days = np.arange(50, np.nanmax(tofs/days2sec), 50)
dv_levels = np.arange(2, cutoff_dv + 0.5, 0.5)

# First plot: C₃, v∞ and TOF   
fig1, ax1 = plt.subplots(figsize=(10, 10))

# Contour  for C3 
c0 = ax1.contourf(X, Y, C3, levels=c3_levels, cmap='magma')
fig1.colorbar(c0, ax=ax1, label=r'C3 ($\mathrm{km^2/s^2}$)')

# Contour for v∞
c1 = ax1.contour(X, Y, v_infs, levels=v_inf_levels, colors='deepskyblue', linewidths=1.5)
ax1.clabel(c1, fmt='%i')

# Contour for TOF
c2 = ax1.contour(
    X, Y, tofs/days2sec,
    levels=tof_levels_days,
    colors='green',
    linestyles='--',   
    linewidths=1
)
ax1.clabel(c2, fmt='%i', colors='green')

# Legend 
ax1.plot([], [], color='deepskyblue', label=r'$V_{\infty}$ (km/s)')
ax1.plot([], [], color='green', linestyle='--', label='TOF (giorni)')

ax1.legend(loc='upper right')
ax1.set_title('Earth–Mars Porkchop Plot (C3, v∞, TOF)')
ax1.set_xlabel('Launch Date')
ax1.set_ylabel('Arrival Date')
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
ax1.yaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
fig1.autofmt_xdate()

plt.show()

# Second plot: total Δv and TOF  

fig2, ax2 = plt.subplots(figsize=(10, 10))

# Contour  for Δv 
cf = ax2.contourf(X, Y, Dv, levels=dv_levels, cmap="magma", alpha=0.85)
cbar = fig2.colorbar(cf, ax=ax2)
cbar.set_label("Δv totale (km/s)")
cbar.set_ticks(dv_levels) 

# Contour for TOF
c_tof = ax2.contour(X, Y, tofs/days2sec, levels=tof_levels_days, colors='green', linestyles='--', linewidths=1)
ax2.clabel(c_tof, fmt='%i', colors='green')

# Legend
ax2.plot([], [], color='green', linestyle='--', label='TOF (giorni)')
ax2.legend(loc='upper right')

# Title and axis
ax2.set_title('Earth–Mars Porkchop Plot, Δv totale')
ax2.set_xlabel('Launch Date')
ax2.set_ylabel('Arrival Date')
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
ax2.yaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
fig2.autofmt_xdate()

# Light grid for better readability
ax2.grid(True, alpha=0.3)

plt.show()

### PLOT THE MINIMUM DV TRAJECTORY

# Find the minimum Δv and the corresponding launch and arrival dates
min_dv_idx = np.unravel_index(np.nanargmin(Dv), Dv.shape)  # 
min_dv = Dv[min_dv_idx]  
min_launch_et = launch_ets[min_dv_idx[1]]  
min_arrival_et = arrival_ets[min_dv_idx[0]] 
min_tof = min_arrival_et - min_launch_et

# Conversion ET -> Date UTC
min_launch_date = datetime.fromisoformat(spice.et2utc(min_launch_et, 'ISOC', 0))
min_arrival_date = datetime.fromisoformat(spice.et2utc(min_arrival_et, 'ISOC', 0))

print(f"Minimum Delta V: {min_dv:.2f} km/s")
print(f"Launch Date: {min_launch_date}")
print(f"Arrival Date: {min_arrival_date}")
print(f"TOF: {min_tof/days2sec:.2f} days")

# Solve the Lambert problem corresponding to the minimum Δv
state_earth, _ = spice.spkezr('EARTH BARYCENTER', min_launch_et, 'ECLIPJ2000', 'NONE', 'SUN')
r_earth = state_earth[:3] 
state_mars, _ = spice.spkezr('MARS BARYCENTER', min_arrival_et, 'ECLIPJ2000', 'NONE', 'SUN')
r_mars = state_mars[:3]  


lambert = pk.lambert_problem(r_earth, r_mars, min_tof, mu_sun)
v1 = np.array(lambert.get_v1()[0]) 
v2 = np.array(lambert.get_v2()[0]) 

# Propagate trajectory
n_points = 100  # Number of points to plot along the trajectory
times = np.linspace(0, min_tof, n_points)  # Time steps from departure to arrival
traj_positions = np.zeros((n_points, 3))  # Vector containing the trajectory positions

# Orbit propagation
for i, t in enumerate(times):
    r, _ = propagate_lagrangian(r_earth, v1, t, mu_sun)
    traj_positions[i] = r

# Extract the Earth and Mars orbits
earth_positions = np.zeros((n_points, 3))
mars_positions = np.zeros((n_points, 3))
for i, t in enumerate(np.linspace(min_launch_et, min_arrival_et, n_points)):
    state_earth, _ = spice.spkezr('EARTH BARYCENTER', t, 'ECLIPJ2000', 'NONE', 'SUN')
    state_mars, _ = spice.spkezr('MARS BARYCENTER', t, 'ECLIPJ2000', 'NONE', 'SUN')
    earth_positions[i] = state_earth[:3]
    mars_positions[i] = state_mars[:3]

# Transfer orbit plot
fig, ax = plt.subplots(figsize=(10, 8))
ax.plot(earth_positions[:, 0], earth_positions[:, 1], 'b-', label='Earth Orbit')
ax.plot(mars_positions[:, 0], mars_positions[:, 1], 'r-', label='Mars Orbit')

ax.plot(traj_positions[:, 0], traj_positions[:, 1], 'g--', label='Transfer Trajectory')

ax.plot(0, 0, 'y*', markersize=15, label='Sun')  
ax.plot(r_earth[0], r_earth[1], 'bo', label='Earth at Launch')
ax.plot(r_mars[0], r_mars[1], 'ro', label='Mars at Arrival')

ax.set_xlabel('X (km, ECLIPJ2000)')
ax.set_ylabel('Y (km, ECLIPJ2000)')
ax.set_title(f'Earth–Mars Transfer Trajectory (Min Δv = {min_dv:.2f} km/s)')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_aspect('equal')  

plt.show()
