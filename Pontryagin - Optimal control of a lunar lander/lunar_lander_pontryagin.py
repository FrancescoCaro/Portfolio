# Fuel-optimal lunar landing (maximum final mass) via Pontryagin's Maximum Principle.
# 1-D vertical motion, constant gravity, bounded thrust 0 <= u <= 1.
#
# State    : x1 = altitude h of the centre of mass, x2 = velocity v, x3 = mass m
# Control  : u in [0, 1]   (T = Tmax*u,  mdot = -Tmax*u/c,  c = g0*Isp)
#
# Cost     : maximise the final mass  J = x3(tf)
#
# Hamiltonian (to be maximised over u):
#     H = p1*x2 + p2*(Tmax*u/x3 - g) - p3*Tmax*u/c
# Costate equations:
#     p1' = -dH/dx1 = 0
#     p2' = -dH/dx2 = -p1
#     p3' = -dH/dx3 = p2*Tmax*u/x3^2
# Switching function (H = H1(u) + H2, with H1 = Tmax*u*theta):
#     theta = p2/x3 - p3/c   ->   u = 1 if theta > 0,  u = 0 if theta < 0
#     (theta = 0 over a finite interval would be a singular arc; it is shown
#      analytically to be impossible here: it would give theta = -1/c != 0)
#
# Boundary / transversality conditions (tf free, final mass free):
#     x1(tf) = dx1   (height of the CoM above the landing legs)
#     x2(tf) = 0     (soft landing)
#     p3(tf) = 1
#     H(tf)  = 0
#
# Shooting scheme:
#     [p1(0), p2(0), p3(0), tf]  ->  ODE integration  ->  residual vector
#     [x1f - dx1, x2f, p3f - 1, Hf]  ->  root finder (least squares)  ->  repeat until J = 0
#
# The ODE is integrated with event detection on theta = 0, so the bang-bang
# switch is located exactly and the residuals stay smooth w.r.t. the unknowns.
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, brentq

# ------------------------------------------------------------------ parameters
g = 1.62             # lunar gravity [m/s^2]
g0 = 9.80665         # standard gravity [m/s^2]
Isp = 311.0          # specific impulse [s]
c = g0 * Isp         # effective exhaust velocity [m/s]
Tmax = 6000.0        # maximum thrust [N]

# initial conditions (altitude of the CoM, vertical velocity, mass)
h0, v0, m0 = 1500.0, -40.0, 2000.0   # [m], [m/s], [kg]
dx1 = 1.5                            # final CoM altitude = CoM height above legs [m]

# scaling factors used to non-dimensionalise the residuals
SCALE = dict(h=100.0, v=10.0, H=10.0)


# -------------------------------------------------------------------- dynamics
def theta_fun(y):
    # Switching function theta = p2/x3 - p3/c.
    _, _, m, _, p2, p3 = y
    return p2 / m - p3 / c


def rhs(t, y, u):
    # Right-hand side of the state + costate system for a given throttle u.
    h, v, m, p1, p2, p3 = y
    return [
        v,                       # x1' = x2
        Tmax * u / m - g,        # x2' = Tmax*u/x3 - g
        -Tmax * u / c,           # x3' = -Tmax*u/c
        0.0,                     # p1' = 0
        -p1,                     # p2' = -p1
        p2 * Tmax * u / m**2,    # p3' = p2*Tmax*u/x3^2
    ]


def hamiltonian(y, u):
    h, v, m, p1, p2, p3 = y
    return p1 * v + p2 * (Tmax * u / m - g) - p3 * Tmax * u / c


def propagate(z, n_pts=400):
    # Integrate state and costates from t=0 to tf with the bang-bang law
    # u = 1 if theta > 0 else 0.  z = [p1(0), p2(0), p3(0), tf].
    # Returns time, solution array (6 x N), throttle history u(t), final u.
    p1_0, p2_0, p3_0, tf = z
    y = np.array([h0, v0, m0, p1_0, p2_0, p3_0], dtype=float)
    t = 0.0
    T_all, Y_all, U_all = [], [], []
    u = 1.0 if theta_fun(y) > 0 else 0.0
    n_switch = 0

    while t < tf - 1e-12 and n_switch < 6:
        # terminal event: theta crosses zero (only in the direction that
        # actually flips the current throttle level)
        def ev(tt, yy, u_arg):
            return theta_fun(yy)
        ev.terminal = True
        ev.direction = -1 if u == 1.0 else 1

        sol = solve_ivp(rhs, (t, tf), y, args=(u,), events=ev, method="DOP853",
                        rtol=1e-11, atol=1e-11, dense_output=True)
        tt = np.linspace(t, sol.t[-1], max(n_pts // 4, 5))
        T_all.append(tt)
        Y_all.append(sol.sol(tt))
        U_all.append(np.full_like(tt, u))
        t, y = sol.t[-1], sol.y[:, -1].copy()
        if sol.status == 1:          # event reached -> switch throttle level
            u = 1.0 - u
            n_switch += 1
        else:
            break

    return np.concatenate(T_all), np.hstack(Y_all), np.concatenate(U_all), u


def residuals(z):
    # Boundary-condition residuals (the 'J' of the block diagram).
    if z[3] <= 1.0:
        return np.full(4, 1e3)
    try:
        t, Y, U, u_end = propagate(z)
    except Exception:
        return np.full(4, 1e3)
    yf = Y[:, -1]
    if yf[2] <= 50.0 or t[-1] < z[3] - 1e-6:   # unphysical mass / failed integration
        return np.full(4, 1e3)
    return np.array([
        (yf[0] - dx1) / SCALE["h"],            # x1(tf) = dx1
        yf[1] / SCALE["v"],                    # x2(tf) = 0
        yf[5] - 1.0,                           # p3(tf) = 1
        hamiltonian(yf, u_end) / SCALE["H"],   # H(tf) = 0
    ])


# --------------------------------------------------------------- initial guess
def tf_guess():
    # Estimate of tf from a 'suicide burn' trajectory: free-fall coasting
    # followed by a constant-deceleration full-thrust burn (mass held at m0).
    a = Tmax / m0 - g
    # coast until the stopping distance v_s^2/(2a) equals the remaining height
    f = lambda ts: (h0 - dx1 + v0 * ts - 0.5 * g * ts**2) - (v0 - g * ts) ** 2 / (2 * a)
    ts = brentq(f, 0, 500)
    return ts + abs(v0 - g * ts) / a


def solve_shooting(verbose=True, n_tries=200, seed=0):
    # Solve the two-point boundary value problem by shooting + multistart.
    tf_g = tf_guess()
    rng = np.random.default_rng(seed)
    best = None
    for k in range(n_tries):
        if k == 0:
            z0 = np.array([-1e-3, 1.0, 1.0, tf_g])       # deterministic first guess
        else:                                            # random restarts
            z0 = np.array([
                rng.uniform(-5e-2, 5e-2),
                rng.uniform(-1, 1) * 10 ** rng.uniform(0, 3),
                rng.uniform(0.5, 1.5),
                tf_g * rng.uniform(0.8, 1.3),
            ])
        try:
            # Levenberg-Marquardt least squares on the 4 residuals (4 unknowns)
            sol = least_squares(residuals, z0, method="lm", xtol=1e-14, ftol=1e-14,
                                gtol=1e-14, max_nfev=400)
        except Exception:
            continue
        cost, x = sol.cost, sol.x
        if best is None or cost < best[0]:
            best = (cost, x)
        if cost < 1e-16:
            if verbose:
                print(f"Converged at attempt {k}, cost = {cost:.2e}")
            break
    if verbose:
        print("Final residuals:", residuals(best[1]))
    return best[1]


def save_figure(fig, name, dpi=140):
    # Save a figure next to this script (falls back to the home directory if
    # that folder is read-only). Never crashes the run if saving is impossible.
    folders = [os.path.dirname(os.path.abspath(__file__)), os.path.expanduser("~")]
    for folder in folders:
        path = os.path.join(folder, name)
        try:
            fig.savefig(path, dpi=dpi)
            print(f"Saved figure: {path}")
            return
        except OSError:
            continue
    print(f"Warning: could not save {name} (no writable folder found)")


# ------------------------------------------------------------------------ main
if __name__ == "__main__":
    z = solve_shooting()
    p1_0, p2_0, p3_0, tf = z
    print(f"\np1(0) = {p1_0:.6e}\np2(0) = {p2_0:.6e}\np3(0) = {p3_0:.6e}\ntf    = {tf:.4f} s")

    t, Y, U, _ = propagate(z, n_pts=2000)
    h, v, m, p1, p2, p3 = Y
    theta = p2 / m - p3 / c
    H = p1 * v + p2 * (Tmax * U / m - g) - p3 * Tmax * U / c

    # switching instants
    idx = np.where(np.diff(U) != 0)[0]
    t_sw = t[idx + 1] if len(idx) else []
    print("Switching times [s]:", t_sw)
    print(f"Final mass = {m[-1]:.3f} kg  (propellant used = {m0 - m[-1]:.3f} kg)")
    print(f"h_f = {h[-1]:.6f} m (target {dx1}),  v_f = {v[-1]:.2e} m/s,  "
          f"p3_f = {p3[-1]:.6f},  H_f = {H[-1]:.2e}")

    def mark_switches(axes):
        for a_ in np.atleast_1d(axes).flatten():
            a_.grid(alpha=0.3)
            for i, ts in enumerate(t_sw):
                a_.axvline(ts, color="r", ls="--", lw=0.8,
                           label="Throttle switch" if i == 0 else None)

    # ---- figure 1: states and costates
    fig, ax = plt.subplots(3, 2, figsize=(12, 10), sharex=True)
    panels = [("Altitude $x_1$ [m]", h), ("Velocity $x_2$ [m/s]", v), ("Mass $x_3$ [kg]", m),
              ("Costate $p_1$ [-]", p1), ("Costate $p_2$ [-]", p2), ("Costate $p_3$ [-]", p3)]
    for a_, (name, sig) in zip(ax.flatten(order="F"), panels):
        a_.plot(t, sig, lw=2)
        a_.set_ylabel(name)
    mark_switches(ax)
    ax[0, 0].legend()
    ax[2, 0].set_xlabel("Time [s]")
    ax[2, 1].set_xlabel("Time [s]")
    ax[0, 0].set_title("State variables")
    ax[0, 1].set_title("Costates (Lagrange multipliers)")
    fig.suptitle("Optimal lunar landing - Pontryagin shooting solution")
    fig.tight_layout()
    save_figure(fig, "state_costates.png")

    # ---- figure 2: control, switching function, Hamiltonian
    fig2, ax2 = plt.subplots(3, 1, figsize=(8, 7), sharex=True)
    ax2[0].step(t, U, where="post", lw=2)
    ax2[0].set_ylabel("Throttle $u$ [-]")
    ax2[0].set_ylim(-0.1, 1.1)
    ax2[0].set_title("Optimal control, switching function and Hamiltonian")
    ax2[1].plot(t, theta, lw=2)
    ax2[1].axhline(0, color="k", lw=0.8)
    ax2[1].set_ylabel(r"Switching fn $\theta$")
    ax2[2].plot(t, H, lw=2)
    ax2[2].set_ylabel("Hamiltonian $H$")
    ax2[2].set_xlabel("Time [s]")
    mark_switches(ax2)
    ax2[0].legend()
    fig2.tight_layout()
    save_figure(fig2, "control_switching_H.png")
    plt.show()
