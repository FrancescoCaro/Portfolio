# Fuel-Optimal Lunar Landing via Pontryagin's Maximum Principle

![Python](https://img.shields.io/badge/python-3.9%2B-blue)
![NumPy](https://img.shields.io/badge/NumPy-%E2%9C%94-013243)
![SciPy](https://img.shields.io/badge/SciPy-%E2%9C%94-8CAAE6)
![Matplotlib](https://img.shields.io/badge/Matplotlib-%E2%9C%94-11557c)

Indirect optimal-control solution of a **vertical lunar soft landing that maximises the final mass** (i.e. minimises propellant), obtained by applying **Pontryagin's Maximum Principle** and solving the resulting two-point boundary value problem with a **shooting method** in Python.

The solver recovers the classic **bang-bang** structure, *coast, then full-thrust burn*, together with the full history of states, costates (Lagrange multipliers), switching function and Hamiltonian.

![State and costates](state_costates.png)

## Highlights

- Indirect method: the optimality conditions are derived analytically, and the numerical solver only has to find the unknown initial costates and final time.
- Analytical proof that **singular arcs cannot occur**, so the optimal control is purely bang-bang.
- Exact handling of the throttle discontinuity through **ODE event detection** on the switching function, which keeps the shooting residuals smooth.
- Robust **multistart** strategy for the notoriously sensitive costate initial guess.
- Result cross-checked against an independent analytical relation (see [Verification](#verification)).

## Problem statement

**State**

| Variable | Meaning |
|---|---|
| $x_1 = h$ | altitude of the centre of mass |
| $x_2 = v$ | vertical velocity |
| $x_3 = m$ | mass |

**Control:** throttle $u \in [0, 1]$, with thrust $T = T_{max}\,u$ and exhaust velocity $c = g_0 I_{sp}$.

**Dynamics** (flat surface, constant gravity $g$, motion along the vertical only)

$$
\dot x_1 = x_2, \qquad
\dot x_2 = \frac{T_{max}\,u}{x_3} - g, \qquad
\dot x_3 = -\frac{T_{max}\,u}{c}
$$

**Objective:** maximise the final mass, $J = x_3(t_f)$.

**Boundary conditions**

- $x(t_0) = (h_0, v_0, m_0)$ known.
- $t_f$ **free**, final mass free.
- $x_1(t_f) = \Delta x_1$, the height of the centre of mass above the landing legs.
- $x_2(t_f) = 0$ (soft landing).

## Optimal control theory

Hamiltonian, to be maximised over $u$:

$$
H = p_1 x_2 + p_2\left(\frac{T_{max}\,u}{x_3} - g\right) - p_3 \frac{T_{max}\,u}{c}
$$

Costate equations:

$$
\dot p_1 = 0, \qquad \dot p_2 = -p_1, \qquad \dot p_3 = p_2 \frac{T_{max}\,u}{x_3^2}
$$

Transversality conditions (two free parameters, $t_f$ and $x_3(t_f)$):

$$
p_3(t_f) = 1, \qquad H(t_f) = 0
$$

Writing $H = H_1(u) + H_2$ with $H_1 = T_{max}\,u\,\theta$, the **switching function** is

$$
\theta = \frac{p_2}{x_3} - \frac{p_3}{c}
\qquad\Rightarrow\qquad
u = \begin{cases} 1 & \theta > 0 \\ 0 & \theta < 0 \end{cases}
$$

**No singular arcs.** Suppose $\theta \equiv 0$ on a finite interval. Then $\dot\theta = -p_1/x_3 = 0$, so $p_1 = 0$ and $p_2$ is constant. Since $H \equiv 0$, this forces $p_2 = 0$ and hence $p_3 = 1$. Substituting back gives $\theta = -1/c \neq 0$, a contradiction. The control is therefore bang-bang, and for a landing a single switch from $u = 0$ to $u = 1$ is expected.

## Numerical method

Unknowns: $z = [\,p_1(0),\ p_2(0),\ p_3(0),\ t_f\,]$.

```
z --> integrate state + costates (ODE) --> final-time residuals --> J = 0 ?
^                                                                    |
|------------------- least-squares update <------------------- no ---|
```

The four final conditions are stacked into a residual vector:

$$
\big[\,x_1(t_f) - \Delta x_1,\;\; x_2(t_f),\;\; p_3(t_f) - 1,\;\; H(t_f)\,\big] = 0
$$

Implementation notes:

- **Integration:** `scipy.integrate.solve_ivp` with DOP853 and tight tolerances ($10^{-11}$).
- **Discontinuous control:** each integration segment runs with a constant $u$ and stops on a terminal event at $\theta = 0$. The throttle is then flipped and the integration restarts from the exact switching state. This locates the switch precisely and keeps the residuals smooth with respect to the unknowns, which the root finder needs.
- **Solver:** `scipy.optimize.least_squares` (Levenberg-Marquardt), with the residuals non-dimensionalised through scale factors.
- **Initial guess:** $t_f$ comes from an analytical "suicide burn" estimate (free-fall coast followed by a constant-deceleration burn). The costates have no such physical estimate, so a **random multistart** explores different orders of magnitude until a run converges.
- **Robustness:** unphysical trial points (negative time, mass collapse, failed integration) return a large penalty instead of crashing the solver.

## Results

Default scenario:

| Parameter | Value |
|---|---|
| Lunar gravity $g$ | 1.62 m/s² |
| Specific impulse $I_{sp}$ | 311 s |
| Maximum thrust $T_{max}$ | 6000 N |
| Initial altitude $h_0$ | 1500 m |
| Initial velocity $v_0$ | -40 m/s |
| Initial mass $m_0$ | 2000 kg |
| Final CoM height $\Delta x_1$ | 1.5 m |

| Result | Value |
|---|---|
| Landing time $t_f$ | 47.38 s |
| Throttle switch (coast → full thrust) | $t \approx 9.20$ s |
| Propellant used | 75.11 kg |
| Final mass | 1924.89 kg |
| Initial costates $(p_1, p_2, p_3)$ | $(-1.823 \cdot 10^{-2},\ 0.4502,\ 0.9423)$ |

All final conditions are met to numerical precision: residuals around $10^{-14}$, $p_3(t_f) = 1$, $H(t_f) \approx 0$.

How to read the plots:

- **Costates:** $p_1$ is constant, $p_2$ is linear (since $\dot p_2 = -p_1$), and $p_3$ is constant during the coast and grows to 1 during the burn.
- **Control:** $\theta$ crosses zero once, so the throttle switches once from 0 to 1.
- **Hamiltonian:** it stays at the $10^{-14}$ level for the whole trajectory (the second plot below), as required by the free final time condition. The visible fluctuations are just floating-point noise.

![Control, switching function and Hamiltonian](control_switching_H.png)

## Verification

An independent analytical check on the propellant: after the switch the engine runs at full thrust, so the propellant used is

$$
\Delta m = \frac{T_{max}\,(t_f - t_{switch})}{c} = \frac{6000 \cdot 38.18}{9.80665 \cdot 311} \approx 75.11\ \text{kg},
$$

which matches the mass difference obtained from the integrated trajectory.

## Getting started

Requirements: Python 3.9+.

```bash
git clone <your-repository-url>
cd <your-repository-folder>
pip install -r requirements.txt
python lunar_lander_pontryagin.py
```

The script prints the solution and saves `state_costates.png` and `control_switching_H.png` next to the script (falling back to the home directory if that folder is read-only).

All physical parameters and initial conditions are defined at the top of `lunar_lander_pontryagin.py`, so other landers or scenarios can be tried by editing a few lines.

### Convergence tips

Shooting methods are sensitive to the initial costate guess. If the solver does not converge after changing parameters:

- increase `n_tries` in `solve_shooting` for more random restarts;
- widen the random search ranges for $p_1(0)$, $p_2(0)$, $p_3(0)$;
- make sure the problem is feasible: $T_{max}/m_0 > g$ so that the lander can decelerate, and enough altitude is available to stop.

## Assumptions and limitations

- Vertical 1-D motion, constant gravity, no atmosphere.
- No minimum throttle and no dry-mass constraint.
- The solution structure is assumed to be bang-bang (the code allows up to 6 switches).
- The solution is an extremal satisfying the necessary conditions of the Maximum Principle. For this problem it coincides with the known optimal structure.

## Possible extensions

- Minimum throttle level and dry-mass constraint.
- Altitude-dependent gravity or a 2-D descent with a pitch-angle control.
- Time-optimal or mixed cost functions.
- A homotopy or direct-collocation method to generate better initial guesses for the costates.

## Repository structure

| File | Description |
|---|---|
| `lunar_lander_pontryagin.py` | Solver and plotting script |
| `requirements.txt` | Python dependencies |
| `state_costates.png` | State and costate histories |
| `control_switching_H.png` | Throttle, switching function and Hamiltonian |

## References

- L. S. Pontryagin, V. G. Boltyanskii, R. V. Gamkrelidze, E. F. Mishchenko, *The Mathematical Theory of Optimal Processes*, 1962.
- J. S. Meditch, "On the problem of optimal thrusting programs for a lunar soft landing", *IEEE Transactions on Automatic Control*, 1964.

## Author

**Your Name** · [LinkedIn](https://www.linkedin.com/in/your-profile) · [GitHub](https://github.com/your-username)

## License

Released under the MIT License. Add a `LICENSE` file to the repository.
