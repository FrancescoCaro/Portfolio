
Earth-Mars Porkchop Plot Generator
This repository contains a single Python script (Porkchopplot_generator_thesis.py) that generates porkchop plots for ballistic Earth–Mars transfers and identifies the minimum-Δv trajectory within user-defined launch and arrival windows.

The code was developed as part of a master's thesis in aerospace engineering and combines SPICE (via spiceypy) for high-precision planetary ephemerides with pykep for solving Lambert's problem.

What the Code Does
The script performs a brute-force scan over a grid of possible launch and arrival dates and, for each valid pair:

Retrieves heliocentric position and velocity of Earth (at departure) and Mars (at arrival) using NASA's SPICE kernels.
Solves Lambert's problem (pykep) to find the transfer orbit that connects the two positions in the given time-of-flight (TOF).
Computes:
C3 (launch characteristic energy) = ‖v₁ − v⊕‖²
v∞ at Mars arrival = ‖v₂ − v♂‖
Δv at Earth departure = √(C3 + v_escape²) − v_circ,LEO
Δv at Mars arrival (for circularization into 200 km LMO)
Total Δv = Δv_departure + Δv_arrival
Applies user-defined cutoffs (C3 ≤ 50 km²/s², v∞ ≤ 12 km/s, total Δv ≤ 20 km/s) – values above cutoff are clipped for plotting purposes.
Stores results in matrices and generates three figures:
Porkchop Plot 1: Filled contours of C3, blue contours of v∞ (km/s), green dashed contours of TOF (days)
Porkchop Plot 2: Filled contours of total Δv (km/s) with green dashed TOF contours – this is the classic "mission design" porkchop plot
Optimal Transfer Trajectory: Heliocentric plot (ECLIPJ2000) showing Earth orbit, Mars orbit, and the minimum-Δv transfer arc during the flight period, with Sun, departure Earth, and arrival Mars marked.
At the end it prints the minimum total Δv found and the corresponding launch date, arrival date, and TOF.

Default Parameters (2020–2021 Window)
Launch window: 2020-01-01 → 2020-12-31
Arrival window: 2020-10-01 → 2021-12-12
Time step: 1 day
Parking orbit: LEO at 300 km altitude
Mars capture orbit: 200 km circular
Cutoffs: C3 = 50 km²/s², v∞ = 12 km/s, Δv = 20 km/s
These can be easily changed by editing the variables near the top of the script.

Requirements
Bash
pip install numpy matplotlib spiceypy pykep
You also need three NAIF/SPICE kernels (paths are currently hardcoded):

de440s.bsp        → small-body-extended planetary ephemeris (includes Earth & Mars barycenters 1950–2050)
naif0012.tls       → leap seconds kernel
Gravity.tpc        → textual PCK with GM values (or any recent pck file containing BODY10_GM, BODY399_GM, BODY499_GM)
Change the spice.furnsh() paths to match your local setup.

How to Run
Bash
python Porkchopplot_generator_thesis.py
The script will:

Compute the full grid (≈130 000 Lambert problems → takes ~30–60 seconds on a modern laptop)
Display the three matplotlib figures
Print the optimal solution to the console
Example output (for the default 2020–2021 window):

text
Minimum Delta V: 5.87 km/s
Launch Date: 2020-08-25 00:00:00
Arrival Date: 2021-03-11 00:00:00
TOF: 198.00 days
Customizing the Script
To analyze a different launch opportunity (e.g. g. the 2033–2035 synodic period), simply modify these lines:

Python
launch_start_day = "2033-01-01 00:00:00 UTC"
launch_end_day   = "2035-12-31 00:00:00 UTC"
arrival_start_day = "2033-09-01 00:00:00 UTC"
arrival_end_day   = "2036-06-01 00:00:00 UTC"
You can also change orbit altitudes, cutoff values, or the time step (smaller step → higher resolution but much longer runtime).
