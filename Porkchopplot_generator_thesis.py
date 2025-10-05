import numpy as np
import spiceypy as spice
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import pykep as pk
from pykep import propagate_lagrangian

# Carico i kernel di NASA/NAIF 
spice.furnsh("/mnt/c/Users/Francesco/Desktop/Tesi/Gravity.tpc")  # Contiene i valori di GM
spice.furnsh("/mnt/c/Users/Francesco/Desktop/Tesi/naif0012.tls")  # LSK (Leap Seconds Kernel)
spice.furnsh("/mnt/c/Users/Francesco/Desktop/Tesi/de440s.bsp")    # SPK (Planet Ephemeris Kernel)

# Definiamo C3, v_inf e Dv di cutoff oltre il quale non siamo interessati a effettuare i calcoli
cutoff_c3 = 50
cutoff_dv = 20
cutoff_v_inf =  12

### DEFINIZIONE RANGE GIORNI DI LANCIO E ARRIVO
# Definisco la finestra di lancio
launch_start_day =  "2020-01-01 00:00:00 UTC" 
launch_end_day  = "2020-12-31 00:00:00 UTC"

# Definisco la finestra di arrivo
arrival_start_day = "2020-10-01 00:00:00 UTC"
arrival_end_day  = "2021-12-12 00:00:00 UTC"

step = 1 # Step temporale
days2sec = 86400.0
step_s = step*days2sec

# Conversione estremi in ET (secondi da J2000)
et_launch_start = spice.utc2et(launch_start_day)
et_launch_end   = spice.utc2et(launch_end_day)
et_arrival_start = spice.utc2et(arrival_start_day)
et_arrival_end   = spice.utc2et(arrival_end_day)

# Costruzione dei vettori contenenti le possibili date di lancio e di arrivo
launch_ets = np.arange(et_launch_start, et_launch_end +1, step_s, dtype=float)
arrival_ets = np.arange(et_arrival_start, et_arrival_end +1, step_s, dtype=float)

# Numero di giorni per ogni finestra e combinazioni totali
ls = len(launch_ets)
as_ = len(arrival_ets)
total_c = ls*as_

### COSTRUZIONE ALGORITMO

# Leggi le costanti gravitazionali [km^3]/[s^2]
mu_sun   = spice.gdpool("BODY10_GM", 0, 1)[0]
mu_earth = spice.gdpool("BODY399_GM", 0, 1)[0]
mu_mars  = spice.gdpool("BODY499_GM", 0, 1)[0] 

# Calcolo delle grandezze orbitali rilevanti nel calcolo del Delta V 
r_LEO = (6371 + 300)      # raggio Terra + LEO 300 km [km]
v_esc = np.sqrt(2*mu_earth/r_LEO)        # velocità di fuga per un'orbita circolare di raggio r_LEO [km/s]

r_LMO = (3390 + 200) # raggio dell'orbita di arrivo di Marte [km]
v_esc_M =  np.sqrt(2*mu_mars/r_LMO)  # velocità di fuga sull'orbita circolare attorno a Marte di raggio r_LMO [km/s]


# Inizializzo i vettori contententi le velocità e i TOF
C3 = np.full((as_,ls),np.nan)    # Matrice che conterrà i valori di C3
tofs = np.full((as_,ls),np.nan)   # Matrice che contiene i TOF
Dv = np.full((as_,ls),np.nan)     # Matrice contenete i Dv
v_infs = np.full((as_,ls),np.nan)  # Matrice contenente i v_inf_arr

## Ciclo principale per costruire le matrici TOF, C3, v_inf_arr e DeltaV
for i, launch_date in enumerate(launch_ets):

   # Ottengo la posizione dei pianeti nel sistema di riferimento eliocentrico
    state_earth, lt = spice.spkezr('EARTH BARYCENTER', launch_date, 'ECLIPJ2000', 'NONE', 'SUN') 
    r_earth = state_earth[:3]   # Vettore posizione della Terra nel sdr del Sole [km]
    v_earth = state_earth[3:]   # Vettore velocità della Terra nel sdr del Sole [km/s]

    for j, arrival_date in enumerate(arrival_ets):

        # Salta combinazioni in cui la data di arrivo è precedente alla data di lancio
        if arrival_date <= launch_date:
            continue

        # Calcolo del TOF per ogni combinazione di launch_date e arrival_date in secondi
        TOF = arrival_date - launch_date  
        tofs[j, i] = TOF  # Salvo TOF nella matrice tofs

        state_mars, lt = spice.spkezr('MARS BARYCENTER', arrival_date, 'ECLIPJ2000', 'NONE', 'SUN')
        r_mars = state_mars[:3]     # Vettore posizione di Marte nel sdr del Sole [km]
        v_mars = state_mars[3:]     # Vettore velocità di Marte nel sdr del Sole [km/s]

        try:
            # Risolvo il problema di Lambert tra la Terra e Marte per il TOF calcolato
            sol = pk.lambert_problem(r_earth, r_mars, TOF, mu_sun)
            # Estraggo le velocità di partenza e arrivo dalla soluzione di Lambert
            v1 = np.array(sol.get_v1()[0])  # velocità di immissione da Terra [km/s]
            v2 = np.array(sol.get_v2()[0])  # velocità al momento dell'arrivo su Marte [km/s]
        except:
            continue

        # Calcolo C3 e v_inf_arr
        c3 = (np.linalg.norm(v1 - v_earth))**2  # [km^2/s^2]
        v_inf_arr = np.linalg.norm(v2 - v_mars)  # [km/s]

        # Calcolo del DeltaV richiesto
        dv1 = np.sqrt(c3 + v_esc**2) - np.sqrt(mu_earth/r_LEO)
        dv2 = np.sqrt(v_inf_arr**2 + v_esc_M**2) - np.sqrt(mu_mars/r_LMO)
        dv = dv1 + dv2  # DeltaV totale [km/s]

        # Applico i cutoff per C3, v_inf_arr e DeltaV
        c3 = min(c3, cutoff_c3)
        v_inf_arr = min(v_inf_arr, cutoff_v_inf)
        dv = min(dv, cutoff_dv)

        # Salvo i valori nelle rispettive matrici
        C3[j, i] = c3   
        v_infs[j, i] = v_inf_arr  
        Dv[j, i] = dv

# --- Conversione ET → datetime ---
launch_dates = [datetime.fromisoformat(spice.et2utc(et, 'ISOC', 0)) for et in launch_ets]
arrival_dates = [datetime.fromisoformat(spice.et2utc(et, 'ISOC', 0)) for et in arrival_ets]

### GRAFICO I RISULTATI

# Creazione griglia coordinate 
X, Y = np.meshgrid(launch_dates, arrival_dates)

#  Livelli per i plot 
c3_levels = np.arange(0, cutoff_c3, 5)
v_inf_levels = np.arange(0, cutoff_v_inf, 3)
tof_levels_days = np.arange(50, np.nanmax(tofs/days2sec), 50)
dv_levels = np.arange(2, cutoff_dv + 0.5, 0.5)

# Primo grafico: C3 + v∞ + TOF 
# Primo grafico: C3 (contourf) + v∞ + TOF 
fig1, ax1 = plt.subplots(figsize=(10, 10))

# Contour  per C3 
c0 = ax1.contourf(X, Y, C3, levels=c3_levels, cmap='magma')
fig1.colorbar(c0, ax=ax1, label=r'C3 ($\mathrm{km^2/s^2}$)')

# Contour per v∞
c1 = ax1.contour(X, Y, v_infs, levels=v_inf_levels, colors='deepskyblue', linewidths=1.5)
ax1.clabel(c1, fmt='%i')

# Contour per TOF
c2 = ax1.contour(
    X, Y, tofs/days2sec,
    levels=tof_levels_days,
    colors='green',
    linestyles='--',   # linee tratteggiate
    linewidths=1
)
ax1.clabel(c2, fmt='%i', colors='green')

# Legenda 
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

# Secondo grafico: Δv totale + TOF 

fig2, ax2 = plt.subplots(figsize=(10, 10))

# Contour  per Δv 
cf = ax2.contourf(X, Y, Dv, levels=dv_levels, cmap="magma", alpha=0.85)
cbar = fig2.colorbar(cf, ax=ax2)
cbar.set_label("Δv totale (km/s)")
cbar.set_ticks(dv_levels)  # Mostra i livelli specifici di Δv

# Contour per TOF
c_tof = ax2.contour(X, Y, tofs/days2sec, levels=tof_levels_days, colors='green', linestyles='--', linewidths=1)
ax2.clabel(c_tof, fmt='%i', colors='green')

# Legenda con Δv e TOF
ax2.plot([], [], color='green', linestyle='--', label='TOF (giorni)')
ax2.legend(loc='upper right')

# Titolo e assi
ax2.set_title('Earth–Mars Porkchop Plot, Δv totale')
ax2.set_xlabel('Launch Date')
ax2.set_ylabel('Arrival Date')
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
ax2.yaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
fig2.autofmt_xdate()

# Griglia leggera per leggibilità
ax2.grid(True, alpha=0.3)

plt.show()

### GRAFICO LA TRAIETTORIA CON IL DV MINIMO

# Trovo il deltaV minimo e le rispettive date
min_dv_idx = np.unravel_index(np.nanargmin(Dv), Dv.shape)  # 
min_dv = Dv[min_dv_idx]  
min_launch_et = launch_ets[min_dv_idx[1]]  
min_arrival_et = arrival_ets[min_dv_idx[0]] 
min_tof = min_arrival_et - min_launch_et

# Converto ET -> Data UTC
min_launch_date = datetime.fromisoformat(spice.et2utc(min_launch_et, 'ISOC', 0))
min_arrival_date = datetime.fromisoformat(spice.et2utc(min_arrival_et, 'ISOC', 0))

print(f"Minimum Delta V: {min_dv:.2f} km/s")
print(f"Launch Date: {min_launch_date}")
print(f"Arrival Date: {min_arrival_date}")
print(f"TOF: {min_tof/days2sec:.2f} days")

# Risolvo Lambert corrispondente al minimo Dv
state_earth, _ = spice.spkezr('EARTH BARYCENTER', min_launch_et, 'ECLIPJ2000', 'NONE', 'SUN')
r_earth = state_earth[:3] 
state_mars, _ = spice.spkezr('MARS BARYCENTER', min_arrival_et, 'ECLIPJ2000', 'NONE', 'SUN')
r_mars = state_mars[:3]  


lambert = pk.lambert_problem(r_earth, r_mars, min_tof, mu_sun)
v1 = np.array(lambert.get_v1()[0]) 
v2 = np.array(lambert.get_v2()[0]) 

# Propagate la traiettoria
n_points = 100  # Numero di punti da plottare lungo la traiettoria
times = np.linspace(0, min_tof, n_points)  # Time steps dal lancio all'arrivo
traj_positions = np.zeros((n_points, 3))  # Vettore che contiene le posizioni della traiettoria

# Propagation dell'orbita
for i, t in enumerate(times):
    r, _ = propagate_lagrangian(r_earth, v1, t, mu_sun)
    traj_positions[i] = r

# Estraggo le orbite di Terra e Marte
earth_positions = np.zeros((n_points, 3))
mars_positions = np.zeros((n_points, 3))
for i, t in enumerate(np.linspace(min_launch_et, min_arrival_et, n_points)):
    state_earth, _ = spice.spkezr('EARTH BARYCENTER', t, 'ECLIPJ2000', 'NONE', 'SUN')
    state_mars, _ = spice.spkezr('MARS BARYCENTER', t, 'ECLIPJ2000', 'NONE', 'SUN')
    earth_positions[i] = state_earth[:3]
    mars_positions[i] = state_mars[:3]

# Plot del trasferimento
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
