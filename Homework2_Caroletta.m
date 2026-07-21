clc
clear all
close all
format long e

%% PARAMETERS OF THE SIMULATION

% Orital and planetary parameters
r = 3430; % Radius of the circular orbit around Mercury [km]
GM = 22032; % Mercury gravitational parameter [km^3/s^2]
D = 0.41; % Distance between the spacecraft and the Sun [AU]
Rm = 2440; % Radius of Mercury [km]

% Solar radiation pressure parameters
PHI = 1371; % Solar flux at 1 AU [W/m^2]
c = 299792.458e3; % Speed of light [m/s]
A = 8; % Solar panel area [m^2]
Cs = 0.3; % Solar panel reflectivity coefficient 
b = 2; % Arm of the srp [m]

% Wheels parameters
I = 6500; % Moment of inertia of the spacecraft
Iw = 0.054; % Moment of inertia of the wheels 
M = [-1/2,  1/2,  1/2, -1/2;  -1/2, -1/2,  1/2,  1/2;sqrt(2)/2, sqrt(2)/2, sqrt(2)/2, sqrt(2)/2]; % Mounting matrix
omega_w_max = 4000*2*pi/60; % Maximum angular velocity of the wheels [rad/s]
omega_w0 = [1751;-121;-1000;0]*2*pi/60; % Initial angular velocities of the wheels [rad/s]
Tc_max = 0.511; % Wheel maximum angular momentum [Nm]
Km = 0.118; % [Nm/A]
Kw = 0.008; % [V/rpm]
R = 1.2; % [Ohm]
Pmax = 29; % Wheel maximum power [W]
i_max = Tc_max / Km; % Wheel maximum current [A]

%% SOLAR RADIATION PRESSURE DISTURBANCE TORQUE ESTIMATION

% Vectors in the ICRF 
X_icrf = [1;0;0];   % Position of the Sun
alpha = deg2rad(-120);  % Fixed angle between the solar panels and the Sun
S_icrf = [cos(alpha); sin(alpha); 0]; % Solar panel direction
n_icrf = [-sin(alpha); cos(alpha); 0]; % Direction normal to the solar panels
b_icrf = b * S_icrf; % Vector of the momentum arm

% Calculation of the torque
cos_theta = dot(n_icrf, X_icrf); % Dot product between n_icrf and X_icrf
F_srp_icrf = -(PHI/c)*(1/D)^2 * A * cos_theta * ( (1 - Cs) * X_icrf + 2 * Cs * cos_theta * n_icrf ); % SRP force in the ICRF
T_icrf = cross(b_icrf, F_srp_icrf); % SRP torque in the ICRF
T_dy = T_icrf(3); % Disturbance torque affecting the system [Nm]


%% POINT 1 - SIMULATION FROM INITIAL CONDITION

% Dynamics function
function dX = dynamics(t, X, I, Iw, Kp, Kd, Msrp, theta_g, dtheta_g)
    x1 = X(1); % theta
    x2 = X(2); % theta_dot

    % Error calculation
    err = x1 -  theta_g(t) ;
    derr = x2 - dtheta_g(t) ;

    % Control torque calculation
    Tc_req =- Kp * err - Kd * derr;
    Tc_max = 0.511;
    Tc = max(-Tc_max, min(Tc_max, Tc_req));
    
    dx1 = x2;
    dx2 = (Msrp +Tc) / I;
    dx3 = 0;  
    dx4 = -Tc / Iw; 
    dx5 = Tc / Iw;  
    
    dX = [dx1; dx2; dx3;dx4;dx5];
end

% Functions to define the wanted pointing angle
V = sqrt(GM/r); % Radial velocity of the spacecraft [km/s]
n = V / r; % Orbital mean motion
theta_g = @(t) n.*t; 
dtheta_g = @(t) n .* ones(size(t));


% Parameters of the system
X0 = [deg2rad(0.1);0;omega_w0(1);omega_w0(2);omega_w0(3)]; % Initial conditions
tspan = [0 400]; % Integration period

% Iteration loop to look for the value of Kd and Kp that satisfy the
% requirements
Kp_candidates = logspace(-1, 4, 50);
Kp_opt = NaN;

for k = 1:length(Kp_candidates)
    Kp = Kp_candidates(k);
    Kd = 2*sqrt(Kp*I);

    % Simulation
    options = odeset('RelTol',1.0E-13, 'AbsTol',1.0E-13);
    [t, X] = ode113(@(t, X) dynamics(t, X, I, Iw, Kp, Kd, T_dy ,theta_g,dtheta_g), tspan, X0, options);
    
    % Values extraction
    theta_rad = X(:, 1);
    theta_arcsec = theta_rad * (648000 / pi); 
    theta_dot = X(:, 2);
    omega_rpm_w1 = X(:, 3) * (30/pi); 
    omega_rpm_w2 = X(:, 4) * (30/pi);
    omega_rpm_w3 = X(:, 5) * (30/pi);
    
    % Calculation of the applied torque
    Tc_array = -Kp .* (theta_rad - theta_g(t)) - Kd .* (theta_dot - dtheta_g(t));
    Tc_array = max(-Tc_max, min(Tc_max, Tc_array));
    
    % Distribute the torque to the 3 wheels using the inverse mounting matrix
    Tw1_array = zeros(size(Tc_array));
    Tw2_array = -Tc_array;
    Tw3_array = Tc_array;
    
    % Electrical quantities
    % Wheel 1 
    i_w1 = Tw1_array / Km;
    P_w1 = (Kw .* omega_rpm_w1 + R .* i_w1) .* i_w1;
    
    % Wheel 2 
    i_w2 = Tw2_array / Km;
    P_w2 = (Kw .* omega_rpm_w2 + R .* i_w2) .* i_w2;
    
    % Wheel 3
    i_w3 = Tw3_array / Km;
    P_w3 = (Kw .* omega_rpm_w3 + R .* i_w3) .* i_w3;

    % Check if the wheel constraints are satisfied
    constraint_wheel_P= (max(abs(P_w3))<=Pmax) && (max(abs(P_w2))<=Pmax);
    constraint_wheel_i = (max(abs(i_w3))<=i_max) && (max(abs(i_w2))<=i_max);
    constraint_wheel_torque = (max(abs(Tw2_array))<=Tc_max) && (max(abs(Tw3_array))<=Tc_max);
    constraints_wheel = constraint_wheel_P && constraint_wheel_i &&  constraint_wheel_torque;

    % Check if the system constraints are satisfied 
    err_arcsec = (theta_rad - theta_g(t)) * (648000/pi);
    err_at_300 = interp1(t, err_arcsec, 300);
    idx_ss = t > 300;
    err_ss = err_arcsec(idx_ss);
    constraint_ss = max(abs(err_ss)) <= 4;
    constraint_system = (abs(err_at_300) <= 4) && constraint_ss ;  

    % Check if all the constraints are met
    all_constraints_satisfied = constraints_wheel && constraint_system;
    if all_constraints_satisfied
        Kp_opt = Kp;
        break
        % Exit the loop -> Kp_opt is the value that meets all the
        % constraints
    end
end

Kp = Kp_opt;
Kd = 2*sqrt(Kp*I);


figure('Name', 'Control System Analysis - Constant Disturbance', 'Position', [100, 100, 1400, 700]);


% Top-left: Pointing angle
subplot(2,2,1);
theta_g_arcsec = theta_g(t) * (648000/pi);
err_arcsec     = (theta_rad - theta_g(t)) * (648000/pi);

yyaxis left
plot(t, theta_rad * (648000/pi), 'Color', '#0072BD', 'LineWidth', 1.5, 'DisplayName', '\theta (actual)');
hold on;
plot(t, theta_g_arcsec, '--', 'Color', '#77AC30', 'LineWidth', 1.5, 'DisplayName', '\theta_g (guidance)');
ylabel('Angle [arcsec]');

yyaxis right
plot(t, err_arcsec, 'Color', '#D95319', 'LineWidth', 1.2, 'DisplayName', 'Error');
yline(4,  '--r', '+4 arcsec',  'LabelHorizontalAlignment', 'left','LabelVerticalAlignment', 'top', 'HandleVisibility', 'off');
yline(-4, '--r', '-4 arcsec', 'LabelHorizontalAlignment', 'left','LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
ylabel('Error [arcsec]');
xlabel('Time [s]');
title('Pointing angle: actual vs guidance');
legend('Location', 'best');
grid on; grid minor;

% Top-right: Wheel angular velocities
subplot(2,2,2);
plot(t, omega_rpm_w1, 'Color', '#0072BD', 'LineWidth', 1.5, 'DisplayName', 'Wheel 1');
hold on;
plot(t, omega_rpm_w2, 'Color', '#D95319', 'LineWidth', 1.5, 'DisplayName', 'Wheel 2');
plot(t, omega_rpm_w3, 'Color', '#77AC30', 'LineWidth', 1.5, 'DisplayName', 'Wheel 3');
yline(4000,  '--r', '+4000 rpm', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
yline(-4000, '--r', '-4000 rpm', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xlabel('Time [s]'); ylabel('\omega_w [rpm]');
title('Wheel angular velocities');
legend('Location', 'best');
grid on; grid minor;

% Bottom-left: Wheel currents
subplot(2,2,3);
plot(t, i_w1, 'Color', '#0072BD', 'LineWidth', 1.5, 'DisplayName', 'Wheel 1');
hold on;
plot(t, i_w2, 'Color', '#D95319', 'LineWidth', 1.5, 'DisplayName', 'Wheel 2');
plot(t, i_w3, 'Color', '#77AC30', 'LineWidth', 1.5, 'DisplayName', 'Wheel 3');
yline(i_max,  '--r', sprintf('+%.2f A', i_max), 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
yline(-i_max, '--r', sprintf('-%.2f A', i_max), 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xlabel('Time [s]'); ylabel('i_w [A]');
title('Wheel currents');
legend('Location', 'best');
grid on; grid minor;

% Bottom-right: Wheel powers
subplot(2,2,4);
plot(t, P_w1, 'Color', '#0072BD', 'LineWidth', 1.5, 'DisplayName', 'Wheel 1');
hold on;
plot(t, P_w2, 'Color', '#D95319', 'LineWidth', 1.5, 'DisplayName', 'Wheel 2');
plot(t, P_w3, 'Color', '#77AC30', 'LineWidth', 1.5, 'DisplayName', 'Wheel 3');
yline(Pmax,  '--r', sprintf('+%d W', Pmax), 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
yline(-Pmax, '--r', sprintf('-%d W', Pmax), 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xlabel('Time [s]'); ylabel('P_w [W]');
title('Wheel powers');
legend('Location', 'best');
grid on; grid minor;

%% POINT 2 - DESATURATION TIME

% Estimation of the time between the initial condition and the desaturation
% maneuver in case of no eclipse

% Angular momentums of the wheels at time t(end) = 400 s
Jw2 = Iw*omega_rpm_w2(end)*(pi/30);
Jw3 = Iw*omega_rpm_w3(end)*(pi/30);

% Maximum angular momentum
Jw_sat = Iw * omega_w_max * sign(T_dy);


% Time from the initial time t= 0 and the desaturation maneuver t_des
t_des2 = (Jw_sat - Jw2 )/ T_dy + t(end);
t_des3 = (Jw_sat - Jw3 )/ T_dy + t(end);

fprintf('--- Point 2 --- \n')
fprintf('Desaturation time without eclipse - Wheel 2: t = %.2f s\n', t_des2);
fprintf('Desaturation time without eclipse - Wheel 3: t = %.2f s\n', t_des3);

if t_des2 < t_des3
    fprintf('→ Wheel 2 saturates first at t = %.2f s\n', t_des2);
else
    fprintf('→ Wheel 3 saturates first at t = %.2f s\n', t_des3);
end

% Estimation of eclipse time interval
theta_eclipse_enter = pi - asin(Rm/r); % True anomaly at which the eclipse starts: T_dy=0
theta_eclipse_exit = pi + asin(Rm/r); % True anomaly at which the eclips ends 
t_ecl_enter = theta_eclipse_enter/n;
t_ecl_exit = theta_eclipse_exit/n;
t_ecl = t_ecl_exit - t_ecl_enter;

fprintf('Eclipse entry time: t = %.2f s\n', t_ecl_enter);
fprintf('Eclipse duration:   t = %.2f s\n', t_ecl);

% Check if wheel saturates before first eclipse
if t_ecl_enter >= t_des2
    fprintf('Wheel 2 saturates before first eclipse\n');
    t_sat_wheel2 = t_des2;
else
    fprintf('Wheel 2 does not saturate before first eclipse\n');
    % Orbital parameters
    T_orb    = 2*pi / n;                      % Orbital period [s]
    t_sunlit = T_orb - t_ecl;                 % Sunlit time per orbit [s]

    % Angular momentum accumulated per orbit
    delta_J_per_orbit = T_dy * t_sunlit;

    % Total angular momentum budget remaining
    delta_J_total = Jw_sat - Jw2;

    % First partial sunlit phase: from t(end) to first eclipse entry
    t_first_sunlit = t_ecl_enter - t(end);
    delta_J_first  = T_dy * t_first_sunlit;

    % Remaining budget after first partial sunlit phase
    delta_J_remaining = delta_J_total - delta_J_first;

    % Number of complete orbits needed
    n_orbits_complete = floor(abs(delta_J_remaining) / abs(delta_J_per_orbit));

    % Residual angular momentum after complete orbits
    delta_J_residual = delta_J_remaining - n_orbits_complete * delta_J_per_orbit;

    % Time to saturate in the final partial sunlit phase
    t_final_sunlit = delta_J_residual / T_dy;

    % Total saturation time
    t_sat_wheel2 = t(end) + t_first_sunlit + t_ecl + n_orbits_complete * T_orb + t_final_sunlit;

    fprintf('Wheel 2 saturates after %d complete orbits\n', n_orbits_complete + 1);
    fprintf('Wheel 2 saturation time: t = %.2f s (%.2f orbital periods)\n', t_sat_wheel2, t_sat_wheel2 / T_orb);
end


%% POINT 3 - NADIR POINTING WITH FIXED SOLAR PANELS

alpha_body = deg2rad(-100); % Angle betw

function Tdy = compute_Tdy(theta, PHI, c, D, A, Cs, b, alpha_body)
    % Computes the disturbance torque along Y_body (= Z_ICRF)
    % due to solar radiation pressure.
    %
    % INPUTS:
    %   theta       - current spacecraft attitude angle [rad]
    %   PHI         - solar flux at 1 AU [W/m^2]
    %   c           - speed of light [m/s]
    %   D           - distance from the Sun [AU]
    %   A           - solar panel area [m^2]
    %   Cs          - specular reflectivity coefficient [-]
    %   b           - SRP moment arm [m]
    %   alpha_body  - solar panel offset angle in body frame [rad]
    %
    % OUTPUT:
    %   Tdy         - disturbance torque along Y_body = Z_ICRF [Nm]

    % Sun direction in ICRF
    X_icrf = [1; 0; 0];

    % Solar panel angle in ICRF
    alpha_icrf = theta + alpha_body + deg2rad(90) ;   

    % Solar panel direction and normal in ICRF
    S_icrf = [ cos(alpha_icrf); sin(alpha_icrf); 0];
    n_icrf = [-sin(alpha_icrf); cos(alpha_icrf); 0];

    % Cosine of incidence angle
    cos_theta = dot(n_icrf, X_icrf);

    % SRP force in ICRF
    F_srp = -(PHI/c) * (1/D)^2 * A * cos_theta *((1-Cs)*X_icrf + 2*Cs*cos_theta*n_icrf);

    % Moment arm vector in ICRF
    b_icrf = b * S_icrf;

    % Disturbance torque (Z_ICRF component = Y_body)
    T_icrf = cross(b_icrf, F_srp);
    Tdy = T_icrf(3);
end

% New dynamics function
function dX = dynamics_new(t, X, I, Iw, Kp, Kd, theta_g, dtheta_g, PHI, c, D, A, Cs, b, alpha_body,Rm, r, n)
    x1 = X(1); % theta
    x2 = X(2); % theta_dot

    % Current orbital angle
    theta_orb = mod(n*t, 2*pi);

    % Eclipse check
    theta_enter = pi - asin(Rm/r);
    theta_exit  = pi + asin(Rm/r);

    if theta_orb >= theta_enter && theta_orb <= theta_exit
        Tdy = 0;
    else
        Tdy = compute_Tdy(x1 + deg2rad(90), PHI, c, D, A, Cs, b, alpha_body);
    end

    % Error calculation
    err  = x1 - theta_g(t);
    derr = x2 - dtheta_g(t);

    % Control torque calculation
    Tc_req = -Kp * err - Kd * derr;
    Tc_max = 0.511;
    Tc = max(-Tc_max, min(Tc_max, Tc_req));

    dx1 = x2;
    dx2 = (Tdy + Tc) / I;
    dx3 = 0;
    dx4 = -Tc / Iw;
    dx5 =  Tc / Iw;

    dX = [dx1; dx2; dx3; dx4; dx5];
end


% Parameters of the system
X0 = [deg2rad(0.1);0;omega_w0(1);omega_w0(2);omega_w0(3)]; % Initial conditions
tspan = [0 T_orb/2]; % Integration period

% Iteration loop to look for the value of Kd and Kp that satisfy the
% requirements
Kp_candidates = logspace(-1, 4, 10);
Kp_opt = NaN;

for k = 1:length(Kp_candidates)
    Kp = Kp_candidates(k);
    Kd = 2*sqrt(Kp*I);

    % Simulation
    options = odeset('RelTol',1.0E-13, 'AbsTol',1.0E-13);
    [t, X] = ode113(@(t, X) dynamics_new(t, X, I, Iw, Kp, Kd, theta_g, dtheta_g, ...
        PHI, c, D, A, Cs, b, alpha_body,Rm, r, n), tspan, X0, options);
    
    % Values extraction
    theta_rad = X(:, 1);
    theta_arcsec = theta_rad * (648000 / pi); 
    theta_dot = X(:, 2);
    omega_rpm_w1 = X(:, 3) * (30/pi); 
    omega_rpm_w2 = X(:, 4) * (30/pi);
    omega_rpm_w3 = X(:, 5) * (30/pi);
    
    % Calculation of the applied torque
    Tc_array = -Kp .* (theta_rad - theta_g(t)) - Kd .* (theta_dot - dtheta_g(t));
    Tc_array = max(-Tc_max, min(Tc_max, Tc_array));
    
    % Distribute the torque to the 3 wheels using the inverse mounting matrix
    Tw1_array = zeros(size(Tc_array));
    Tw2_array = -Tc_array;
    Tw3_array = Tc_array;
    
    % Electrical quantities
    % Wheel 1 
    i_w1 = Tw1_array / Km;
    P_w1 = (Kw .* omega_rpm_w1 + R .* i_w1) .* i_w1;
    
    % Wheel 2 
    i_w2 = Tw2_array / Km;
    P_w2 = (Kw .* omega_rpm_w2 + R .* i_w2) .* i_w2;
    
    % Wheel 3
    i_w3 = Tw3_array / Km;
    P_w3 = (Kw .* omega_rpm_w3 + R .* i_w3) .* i_w3;

    % Check if the wheel constraints are satisfied
    constraint_wheel_P= (max(abs(P_w3))<=Pmax) && (max(abs(P_w2))<=Pmax);
    constraint_wheel_i = (max(abs(i_w3))<=i_max) && (max(abs(i_w2))<=i_max);
    constraint_wheel_torque = (max(abs(Tw2_array))<=Tc_max) && (max(abs(Tw3_array))<=Tc_max);
    constraints_wheel = constraint_wheel_P && constraint_wheel_i &&  constraint_wheel_torque;

    % Check if the system constraints are satisfied 
    err_arcsec = (theta_rad - theta_g(t)) * (648000/pi);
    err_at_300 = interp1(t, err_arcsec, 300);
    idx_ss = t > 300;
    err_ss = err_arcsec(idx_ss);
    constraint_ss = max(abs(err_ss)) <= 4;
    constraint_system = (abs(err_at_300) <= 4) && constraint_ss;  

    % Check if all the constraints are met
    all_constraints_satisfied = constraints_wheel && constraint_system;
    if all_constraints_satisfied
        Kp_opt = Kp;
        break
        % Exit the loop -> Kp_opt is the value that meets all the
        % constraints
    end
end

Kp = Kp_opt;
Kd = 2*sqrt(Kp*I);

fprintf('--- Point 3 --- \n')
fprintf('Values after the iteration loop:\n Kp = %2f \n Kd = %2f \n',Kp,Kd)
if isnan(Kp)
    disp('There is no Kp that satisfies all the requirements => A PID control is needed ')
end

%%  PID implementation

% New PID dynamics function
function dX = dynamics_pid(t, X, I, Iw, Kp, Kd, Ki, theta_g, dtheta_g, ...
                            PHI, c, D, A, Cs, b, alpha_body, Rm, r, n)
    x1 = X(1); % theta
    x2 = X(2); % theta_dot
    x6 = X(6); % integral of error

    % Current orbital angle
    theta_orb = mod(n*t, 2*pi);

    % Eclipse check
    theta_enter = pi - asin(Rm/r);
    theta_exit  = pi + asin(Rm/r);

    if theta_orb >= theta_enter && theta_orb <= theta_exit
        Tdy = 0;
    else
        Tdy = compute_Tdy(x1 , PHI, c, D, A, Cs, b, alpha_body);
    end

    % Error calculation
    err  = mod(x1 - theta_g(t) + pi, 2*pi) - pi; 
    derr = x2 - dtheta_g(t);

    % PID control torque
    Tc_req = -Kp * err - Kd * derr - Ki * x6;
    Tc_max = 0.511;
    Tc = max(-Tc_max, min(Tc_max, Tc_req));

    if (Tc_req >= Tc_max && err < 0) || (Tc_req <= -Tc_max && err > 0)
        dx6 = 0;   
    else
        dx6 = err; 
    end



    dx1 = x2;
    dx2 = (Tdy + Tc) / I;
    dx3 = 0;
    dx4 = -Tc / Iw;
    dx5 =  Tc / Iw;

    dX = [dx1; dx2; dx3; dx4; dx5; dx6];
end

% New initial conditions 
X0 = [deg2rad(0.1); 0; omega_w0(1); omega_w0(2); omega_w0(3); 0];

% Iterations on the natural frequency omega_n
wn_candidates = logspace(-3,0 , 100);
Kp_opt = NaN;

for k = 1:length(wn_candidates)

    Kp = 3 * I * wn_candidates(k)^2;
    Kd = 3 * I * wn_candidates(k);
    Ki = I * wn_candidates(k)^3;

    % Simulation
    options = odeset('RelTol',1.0E-13, 'AbsTol',1.0E-13);
    [t, X] = ode113(@(t, X) dynamics_pid(t, X, I, Iw, Kp, Kd, Ki,  theta_g, ...
        dtheta_g, PHI, c, D, A, Cs, b, alpha_body,Rm, r, n), tspan, X0, options);
    
    % Values extraction
    theta_rad = X(:, 1);
    theta_arcsec = theta_rad * (648000 / pi); 
    theta_dot = X(:, 2);
    omega_rpm_w1 = X(:, 3) * (30/pi); 
    omega_rpm_w2 = X(:, 4) * (30/pi);
    omega_rpm_w3 = X(:, 5) * (30/pi);
    
    % Calculation of the applied torque
    integral_err = X(:, 6);
    Tc_array = -Kp .* (theta_rad - theta_g(t)) - Kd .* (theta_dot - dtheta_g(t)) - Ki .* integral_err;
    Tc_array = max(-Tc_max, min(Tc_max, Tc_array));
    
    % Distribute the torque to the 3 wheels using the inverse mounting matrix
    Tw1_array = zeros(size(Tc_array));
    Tw2_array = -Tc_array;
    Tw3_array = Tc_array;
    
    % Electrical quantities
    % Wheel 1 
    i_w1 = Tw1_array / Km;
    P_w1 = (Kw .* omega_rpm_w1 + R .* i_w1) .* i_w1;
    
    % Wheel 2 
    i_w2 = Tw2_array / Km;
    P_w2 = (Kw .* omega_rpm_w2 + R .* i_w2) .* i_w2;
    
    % Wheel 3
    i_w3 = Tw3_array / Km;
    P_w3 = (Kw .* omega_rpm_w3 + R .* i_w3) .* i_w3;

    % Check if the wheel constraints are satisfied
    constraint_wheel_P= (max(abs(P_w3))<=Pmax) && (max(abs(P_w2))<=Pmax);
    constraint_wheel_i = (max(abs(i_w3))<=i_max) && (max(abs(i_w2))<=i_max);
    constraint_wheel_torque = (max(abs(Tw2_array))<=Tc_max) && (max(abs(Tw3_array))<=Tc_max);
    constraints_wheel = constraint_wheel_P && constraint_wheel_i &&  constraint_wheel_torque;

    % Check if the system constraints are satisfied 
    err_arcsec = (theta_rad - theta_g(t)) * (648000/pi);
    err_at_300 = interp1(t, err_arcsec, 300);
    idx_ss = t > 300;
    err_ss = err_arcsec(idx_ss);
    constraint_ss = max(abs(err_ss)) <= 4;
    constraint_system = (abs(err_at_300) <= 4) && constraint_ss;  

    % Check if all the constraints are met
    all_constraints_satisfied = constraints_wheel && constraint_system;
    if all_constraints_satisfied
        Kp_opt = Kp;
        break
        % Exit the loop -> Kp_opt is the value that meets all the
        % constraints
    end
end


fprintf('Values after the iteration loop: \n Kp = %2f \n Kd = %2f \n Ki = %2f \n',Kp,Kd, Ki)

tspan = [0 400];
 % Simulation
    options = odeset('RelTol',1.0E-13, 'AbsTol',1.0E-13);
    [t, X] = ode113(@(t, X) dynamics_pid(t, X, I, Iw, Kp, Kd, Ki,  theta_g, ...
        dtheta_g, PHI, c, D, A, Cs, b, alpha_body,Rm, r, n), tspan, X0, options);
    
    % Values extraction
    theta_rad = X(:, 1);
    theta_arcsec = theta_rad * (648000 / pi); 
    theta_dot = X(:, 2);
    omega_rpm_w1 = X(:, 3) * (30/pi); 
    omega_rpm_w2 = X(:, 4) * (30/pi);
    omega_rpm_w3 = X(:, 5) * (30/pi);
    
    % Calculation of the applied torque
    integral_err = X(:, 6);
    Tc_array = -Kp .* (theta_rad - theta_g(t)) - Kd .* (theta_dot - dtheta_g(t)) - Ki .* integral_err;
    Tc_array = max(-Tc_max, min(Tc_max, Tc_array));
    
    % Distribute the torque to the 3 wheels using the inverse mounting matrix
    Tw1_array = zeros(size(Tc_array));
    Tw2_array = -Tc_array;
    Tw3_array = Tc_array;
    
    % Electrical quantities
    % Wheel 1 
    i_w1 = Tw1_array / Km;
    P_w1 = (Kw .* omega_rpm_w1 + R .* i_w1) .* i_w1;
    
    % Wheel 2 
    i_w2 = Tw2_array / Km;
    P_w2 = (Kw .* omega_rpm_w2 + R .* i_w2) .* i_w2;
    
    % Wheel 3
    i_w3 = Tw3_array / Km;
    P_w3 = (Kw .* omega_rpm_w3 + R .* i_w3) .* i_w3;


figure('Name', 'Control System Analysis - Variable Disturbance', 'Position', [100, 100, 1400, 700]);


% Top-left: Pointing angle
subplot(2,2,1);
theta_g_arcsec = theta_g(t) * (648000/pi);
err_arcsec     = (theta_rad - theta_g(t)) * (648000/pi);

yyaxis left
plot(t, theta_rad * (648000/pi), 'Color', '#0072BD', 'LineWidth', 1.5, 'DisplayName', '\theta (actual)');
hold on;
plot(t, theta_g_arcsec, '--', 'Color', '#77AC30', 'LineWidth', 1.5, 'DisplayName', '\theta_g (guidance)');
ylabel('Angle [arcsec]');

yyaxis right
plot(t, err_arcsec, 'Color', '#D95319', 'LineWidth', 1.2, 'DisplayName', 'Error');
yline(4,  '--r', '+4 arcsec',  'LabelHorizontalAlignment', 'left','LabelVerticalAlignment', 'top', 'HandleVisibility', 'off');
yline(-4, '--r', '-4 arcsec', 'LabelHorizontalAlignment', 'left','LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
ylabel('Error [arcsec]');
xlabel('Time [s]');
title('Pointing angle: actual vs guidance');
legend('Location', 'best');
grid on; grid minor;

% Top-right: Wheel angular velocities
subplot(2,2,2);
plot(t, omega_rpm_w1, 'Color', '#0072BD', 'LineWidth', 1.5, 'DisplayName', 'Wheel 1');
hold on;
plot(t, omega_rpm_w2, 'Color', '#D95319', 'LineWidth', 1.5, 'DisplayName', 'Wheel 2');
plot(t, omega_rpm_w3, 'Color', '#77AC30', 'LineWidth', 1.5, 'DisplayName', 'Wheel 3');
yline(4000,  '--r', '+4000 rpm', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
yline(-4000, '--r', '-4000 rpm', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xlabel('Time [s]'); ylabel('\omega_w [rpm]');
title('Wheel angular velocities');
legend('Location', 'best');
grid on; grid minor;

% Bottom-left: Wheel currents
subplot(2,2,3);
plot(t, i_w1, 'Color', '#0072BD', 'LineWidth', 1.5, 'DisplayName', 'Wheel 1');
hold on;
plot(t, i_w2, 'Color', '#D95319', 'LineWidth', 1.5, 'DisplayName', 'Wheel 2');
plot(t, i_w3, 'Color', '#77AC30', 'LineWidth', 1.5, 'DisplayName', 'Wheel 3');
yline(i_max,  '--r', sprintf('+%.2f A', i_max), 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
yline(-i_max, '--r', sprintf('-%.2f A', i_max), 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xlabel('Time [s]'); ylabel('i_w [A]');
title('Wheel currents');
legend('Location', 'best');
grid on; grid minor;

% Bottom-right: Wheel powers
subplot(2,2,4);
plot(t, P_w1, 'Color', '#0072BD', 'LineWidth', 1.5, 'DisplayName', 'Wheel 1');
hold on;
plot(t, P_w2, 'Color', '#D95319', 'LineWidth', 1.5, 'DisplayName', 'Wheel 2');
plot(t, P_w3, 'Color', '#77AC30', 'LineWidth', 1.5, 'DisplayName', 'Wheel 3');
yline(Pmax,  '--r', sprintf('+%d W', Pmax), 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
yline(-Pmax, '--r', sprintf('-%d W', Pmax), 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xlabel('Time [s]'); ylabel('P_w [W]');
title('Wheel powers');
legend('Location', 'best');
grid on; grid minor;

% Recompute T_dy at each time step using the saved state vector X3

Tdy_array = zeros(size(t));
for k = 1:length(t)
    Tdy_array(k) = compute_Tdy(X(k,1) + deg2rad(90), PHI, c, D, A, Cs, b, alpha_body);
end

% Plot of T_dy over one full orbit 

T_orb = 2*pi/n;                          % Orbital period [s]
t_orbit = linspace(0, T_orb, 10000);     % Time vector over one orbit

% Compute T_dy at each instant (assuming nadir-pointing: theta = n*t)
Tdy_orbit = zeros(size(t_orbit));
for k = 1:length(t_orbit)
    theta_k = n * t_orbit(k);            % nadir-pointing angle at time t
    Tdy_orbit(k) = compute_Tdy(theta_k + deg2rad(90), PHI, c, D, A, Cs, b, alpha_body);
end

fprintf('Mean T_dy over one orbit: %.4e Nm\n', mean(Tdy_orbit));
fprintf('Max  T_dy over one orbit: %.4e Nm\n', max(Tdy_orbit));
fprintf('Min  T_dy over one orbit: %.4e Nm\n', min(Tdy_orbit));

% Plot
figure('Name', 'T_dy over one orbit without considering the eclipse phase', 'Position', [100, 100, 900, 400]);
plot(t_orbit/3600, Tdy_orbit, 'Color', '#0072BD', 'LineWidth', 1.5);
yline(T_dy, '--r', 'T_{dy} point 1', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
yline(0, '--k', 'HandleVisibility', 'off');
yline(mean(Tdy_orbit), '--m', 'Mean T_{dy}', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xlabel('Time [h]');
ylabel('T_{dy} [Nm]');
title('SRP disturbance torque over one full orbit without considering the eclipse');
grid on; grid minor;



% Saturation time estimation
N_orbits = 100;
tspan_long = [0, N_orbits * T_orb];

options_long = odeset('RelTol', 1e-13, 'AbsTol', 1e-13,'Events', @(t,X) saturation_event(t, X, omega_w_max));

function [value, isterminal, direction] = saturation_event(t, X, omega_w_max)
    % Triggers when wheel 2 or wheel 3 reaches saturation
    omega_w2 = abs(X(4));
    omega_w3 = abs(X(5));

    value = min(omega_w_max - omega_w2, omega_w_max - omega_w3);
    isterminal = 1;   % stop integration
    direction  = -1;  % triggered when crossing zero from above
end

[t_long, X_long, t_event, ~, ~] = ode113(@(t,X) dynamics_pid(t, X, I, Iw, Kp, Kd,Ki, ...
    theta_g, dtheta_g, PHI, c, D, A, Cs, b, alpha_body, Rm, r, n), ...
    tspan_long, X0, options_long);

if isempty(t_event)
    fprintf('No saturation occurred within %d orbital periods\n', N_orbits);
else
    fprintf('Saturation time: %.2f s\n', t_event);
    fprintf('Saturation time: %.2f h\n', t_event/3600);
    fprintf('Saturation time: %.2f orbital periods\n', t_event/T_orb);

    % Which wheel saturated
    omega_w2_event = X_long(end, 4) * (30/pi);
    omega_w3_event = X_long(end, 5) * (30/pi);
    fprintf('Wheel 2 speed at saturation: %.2f rpm\n', omega_w2_event);
    fprintf('Wheel 3 speed at saturation: %.2f rpm\n', omega_w3_event);
end