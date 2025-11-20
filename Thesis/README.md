# Earth–Mars Porkchop Plot Generator

**Complete porkchop plot generator** for ballistic Earth → Mars transfers with automatic identification of the **minimum-Δv trajectory**.

Developed as part of a Bachelor's thesis in Aerospace Engineering – La Sapienza (2025)
Author: Francesco  Caroletta

## What the Script Does

The script performs a systematic brute-force search over a user-defined grid of launch and arrival dates. For every valid (launch, arrival) pair it:

1. Retrieves high-precision heliocentric position and velocity of Earth (at departure) and Mars (at arrival) using **NASA SPICE kernels**
2. Solves **Lambert’s problem** with the **pykep** library
3. Computes:
   - **C3** (launch characteristic energy)
   - Arrival **v∞** at Mars
   - Departure **Δv** from a 300 km circular LEO
   - Mars orbit insertion **Δv** into a 200 km circular Low Mars Orbit (LMO)
   - **Total Δv**
4. Applies user-defined cut-offs (C3 ≤ 50 km²/s², v∞ ≤ 12 km/s, total Δv ≤ 20 km/s)
5. Automatically generates **three plots**:
   - Porkchop plot with **C3** (filled contours), **v∞** (blue solid lines) and **TOF** (green dashed lines)
   - Classic porkchop plot with **total Δv** (filled contours) + TOF
   - Heliocentric view of the **optimal transfer trajectory** (minimum Δv) showing Earth orbit, Mars orbit, transfer arc, Sun, departure Earth and arrival Mars

The script also prints to console the **best solution found** (dates, TOF, minimum Δv).

## Default Launch Opportunity (2020–2021 synodic period)

| Parameter                  | Value                              |
|----------------------------|------------------------------------|
| Launch window              | 2020-01-01 → 2020-12-31            |
| Arrival window             | 2020-10-01 → 2021-12-12            |
| Time step                  | **1 day**                          |
| Parking orbit              | LEO 300 km altitude                |
| Mars capture orbit         | Circular 200 km LMO                |
| Typical result             | **≈ 5.87 km/s** total Δv           |
| Example optimal launch     | 2020-08-25                         |
| Example optimal arrival    | 2021-03-11                         |
| Example TOF                | 198 days                           |

## Requirements

```bash
pip install numpy matplotlib spiceypy pykep
```

The script will not run without the following three kernels.
Download and place them in your project folder (or update the paths in the script).
| Kernel    | Direct download link          | Description |
|------------------------|-----------------------------|-------------------|
| `de440s.bsp` | [download here](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp) | Planetary ephemeris (1950–2050)   |
| `naif0012.tls` | [download here](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp) | Leap seconds kernel   |
| `gm_de431.tpc` | [download here](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp) | Gravitational parameters (GM)   |

### Quick one liner to download all three
```bash
wget -c https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp
wget -c https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls
wget -c https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc
```
After downloading, update the three `spice.furnsh()` lines at the top of the script with the correct file paths.

## Customization

To analyze a different Earth–Mars launch opportunity (e.g. the 2033–2035 synodic period), simply modify these variables near the top of the script:

```python
# Launch window
launch_start_day  = "2033-01-01 00:00:00 UTC"
launch_end_day    = "2035-12-31 00:00:00 UTC"

# Arrival window
arrival_start_day = "2033-09-01 00:00:00 UTC"
arrival_end_day   = "2036-06-01 00:00:00 UTC"

# Time resolution (1 = 1 day, 0.5 = 12 hours, etc.)
step = 1  # days
```

Other parameters you can easily change
```python
# Orbits
r_LEO = (6371 + 300)   # Earth departure: 300 km altitude LEO
r_LMO = (3390 + 200)   # Mars arrival: 200 km circular LMO

# Cutoff values (values above these will be clipped in plots)
cutoff_c3     = 50     # km²/s²
cutoff_v_inf  = 12     # km/s
cutoff_dv     = 20     # km/s total Δv
```
