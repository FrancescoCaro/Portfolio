clear all
close all
clc
format long e
set(groot,'defaultLineLineWidth',2)
set(groot,'DefaultAxesFontSize',15)

%% PARAMETERS OF THE SYSTEM

global gm Re A_eff m Cd z0 H_ref theta a Omega nu0 c

c = 299792.458; % Light speed [km/sec]

% Earth constants
gm = 3.986004418e5; % Earth gravitational constant [km^3/sec^2]
Re = 6378;   % Earth radius [km]

% Spacecraft dimensions
A_eff = 25e-6;  % Effective area of the spacecraft [km^2]
m = 1000; % Spacecraft mass [kg]
Cd = 2; % Drag coefficient of the spacecraft [adimensional]

% Constants for drag calculation
z0 = 300;    % [km]
H_ref = 14.7;  % [km]

% Position of the station
theta = 225; % [deg]

% Orbital parameters of the GNSS satellite - E11 satellite of the Galileo
% Constellation
a = 29599.8; % Semi-major axis [km]
e = 0; % Eccentricity -> the orbit is circular
Omega = 77.632; % Right Ascension of the Ascending Node - RAAN [deg]
omega = 0; % Argument of Pericentre [deg]
nu0 = 15.153; % True anomaly at t0 [deg]


%% OBSERVABLES LOADING
obs_file = load('observables_gnss.txt');
epochs = obs_file(:,1);     % Time tags: first column of observable file
range_obs = obs_file(:,2);  % Pseudo-Range observables: second column
rate_obs = obs_file(:,3);   % Pseudo-Range rate observables: third column
N_points = length(epochs);

% Observables plot
figure(1)
subplot(2,1,1)
plot(epochs, range_obs, 'k*-')
title('Pseudo-range observed observables')
xlabel('Time [s]')
ylabel('Range [km]')

subplot(2,1,2)
plot(epochs, rate_obs, 'r*-')
title('Pseudo-range rate observed observables')
xlabel('Time [s]')
ylabel('Range rate [km/s]')

%% DEFINITION OF THE FUNCTIONS OF THE DYNAMICAL SYSTEM

function A = dynamical_matrix(X) 
% Compute the dynamical matrix A starting from the state vector X
% State vector: X = [x, y, vx, vy, rho0, dt, dt_dot]'

    global gm Re A_eff m Cd z0 H_ref 

    % Extract state variables
    x = X(1);
    y = X(2);
    vx = X(3);
    vy = X(4);
    rho0 = X(5);
    dt = X(6);
    dt_dot = X(7);

    % Compute derived quantities
    r = sqrt(x^2 + y^2);
    V = sqrt(vx^2 + vy^2);
    z = r - Re;

    % Atmospheric model
    exp_term = exp(-(z - z0)/H_ref);
    rho = rho0 * exp_term;
    k = Cd * A_eff / (2*m);

    % Initialize 7x7 matrix
    A = zeros(7,7);

    % Row 1
    A(1,3) = 1;

    % Row 2
    A(2,4) = 1;

    drag_pos_coeff = (k * rho * V) / (H_ref * r);

    % Row 3
    A(3,1) = -gm/r^3 + 3*gm*x^2/r^5 + drag_pos_coeff * vx * x;
    A(3,2) =           3*gm*x*y/r^5 + drag_pos_coeff * vx * y;
    A(3,3) = -k * rho * (V + vx^2/V);
    A(3,4) = -k * rho * (vx * vy / V);
    A(3,5) = -k * exp_term * V * vx;

    % Row 4
    A(4,1) =           3*gm*x*y/r^5 + drag_pos_coeff * vy * x;
    A(4,2) = -gm/r^3 + 3*gm*y^2/r^5 + drag_pos_coeff * vy * y;
    A(4,3) = -k * rho * (vx * vy / V);
    A(4,4) = -k * rho * (V + vy^2/V);
    A(4,5) = -k * exp_term * V * vy;

    % Row 6
    A(6,7) = 1;

end

function dX = Model_and_transition(~,X)
    % ODE function to integrate both trajectory and state-transition-matrix
    % State vector: X = [x, y, vx, vy, rho0, dt, dt_dot]'
    
    global gm Re A_eff m Cd z0 H_ref
    
    % Initialization (7 state elements + 49 STM elements)
    dX = zeros(56,1);
    
    % Extract state for readability
    x = X(1);
    y = X(2);
    vx = X(3);
    vy = X(4);
    rho0 = X(5);
    dt = X(6);
    dt_dot = X(7);
    
    % Derived parameters
    r = sqrt(x^2 + y^2);
    V = sqrt(vx^2 + vy^2);
    z = r - Re;
    rho = rho0 * exp(-(z - z0)/H_ref);
    
    % State Integration
    dX(1) = vx;
    dX(2) = vy;
    dX(3) = -gm*x/r^3 - 0.5 * rho * Cd * (A_eff/m) * V * vx;
    dX(4) = -gm*y/r^3 - 0.5 * rho * Cd * (A_eff/m) * V * vy;
    dX(5) = 0;
    dX(6) = dt_dot;
    dX(7) = 0;
    
    % STM Integration
    % STM is a 7x7 matrix -> 49 elements -> indices 8 to 56
    phi = X(8:56);
    PHI = reshape(phi, 7, 7);
    
    A = dynamical_matrix(X);
    dphi = A * PHI;
    
    dX(8:56) = reshape(dphi, 49, 1);
end


function [ x, y, Vx , Vy] = gnss(nu)
    % Compute the velocity component Vx and Vy and positions x and y of the GNSS satellite
    % E11 satellite Galileo from the true anomaly nu in input

    global gm a Omega 
    
    u = Omega + nu;
    vc = sqrt(gm/a);

    % Position 
    x  =  a*cosd(u);
    y  =  a*sind(u);
    % Velocity
    Vx = -vc*sind(u);
    Vy =  vc*cosd(u);

end

function nu = propagate(epochs_vector)
    % Compute the true anomaly nu at each epoch of the input vector epochs
    % vector

    global gm a nu0

    meanMotion = sqrt(gm/a^3); % Mean motion of the orbit [rad/s]
    nu = mod(nu0 + rad2deg(meanMotion * epochs_vector ), 360);
end 

function Obs = CmpObsSet(X,epochs_vector) 
    % Compute the range and range rate observable.
    %
    % Input:
    %   - X: [Nx7] array containing the state vector at the different
    %              epochs
    %
    % Output:
    %   - Obs: array [Nx2] containing pseudo-range and pseudo-range rate
    
    global c    
    
    % initialization
    Obs = zeros(size(X,1),2);
    x = X(:,1);
    y = X(:,2);
    vx = X(:,3);
    vy = X(:,4);
    dt = X(:,6);
    dt_dot = X(:,7);

    % Computation of the position and velocity of the GNSS satellite at
    % each epoch of the observables file
    nu = propagate(epochs_vector);  
    [X_GNSS, Y_GNSS, Vx_GNSS, Vy_GNSS] = gnss(nu);
    
    % computation
    for i=1:size(X,1)
        rho_geo = sqrt((x(i)-X_GNSS(i))^2 + (y(i)-Y_GNSS(i))^2);
        range = rho_geo + c*dt(i);
        rate  = ((x(i)-X_GNSS(i))*(vx(i)-Vx_GNSS(i)) + (y(i)-Y_GNSS(i))*(vy(i)-Vy_GNSS(i))) / (rho_geo + c*dt(i)) + c*dt_dot(i);
        Obs(i,1) = range;
        Obs(i,2) = rate;
    end

end

function H = Htilde(Xt, Xg)
    % Mapping matrix [2x7] for GNSS pseudo-range and pseudo-range rate
    % Input:
    %   Xt - [7x1] or [1x7] spacecraft state [x y vx vy rho0 dt dt_dot]
    %   Xg - [4x1] or [1x4] GNSS state [xg yg vxg vyg]

    global c

    dx  = Xt(1) - Xg(1);
    dy  = Xt(2) - Xg(2);
    dvx = Xt(3) - Xg(3);
    dvy = Xt(4) - Xg(4);

    rho = sqrt(dx^2 + dy^2);
    D   = rho + c*Xt(6);          
    N   = dx*dvx + dy*dvy;

    H = zeros(2,7);

    % Pseudo-range
    H(1,1) = dx/rho;
    H(1,2) = dy/rho;
    H(1,6) = c;

    % Pseudo-range rate
    H(2,1) = dvx/D - N*dx/(rho*D^2);
    H(2,2) = dvy/D - N*dy/(rho*D^2);
    H(2,3) = dx/D;
    H(2,4) = dy/D;
    H(2,6) = -c*N/D^2;
    H(2,7) = c;
end


%% FILTER CODE

% Bias and drift rate initial condition and a priori accuracy
dt0 = 2e-6;
dt_dot0 = 2.1e-8; 
sigdt_ap = 1;
sigdt_dot_ap = 0.001;

% Putative accuracies of the pseudo-range and pseudo-range rate observables
sigma_rho = 1e-4;
sigma_rate = 1e-5;
R = zeros(2,2);
R(1,1) = sigma_rho^2;
R(2,2) = sigma_rate^2;

% Propagated data from the ls-solution of the other 5 state vector
% components at time t0_gnss
X_apriori_5 = [1.563298746857516e+03; -1.045029890058690e+04; 5.264051085239029e+00; 2.439891724000302e+00; 1.708184359913058e-02];
P_apriori_5 = [
    5.625589749769277e-09, 4.641009959688605e-09, 2.910134043432608e-13, 4.820001893280473e-12, 2.732557271636359e-11;
    4.641009959688605e-09, 3.868591054719024e-09, 2.528933571921944e-13, 4.006431295235115e-12, 2.452343372978348e-11;
    2.910134043432608e-13, 2.528933571921944e-13, 2.553744793939443e-17, 2.515403326929971e-16, 1.954263621575061e-15;
    4.820001893280473e-12, 4.006431295235115e-12, 2.515403326929971e-16, 4.169192540834371e-15, 1.665268582479396e-14;
    2.732557271636359e-11, 2.452343372978348e-11, 1.954263621575061e-15, 1.665268582479396e-14, 1.343302653533971e-11;
    ];

% Merge the initial conditions and the covariance matrix
X0 = [X_apriori_5; dt0; dt_dot0 ];
P0 = blkdiag(P_apriori_5, sigdt_ap^2, sigdt_dot_ap^2);

% Adding t0_GNSS to the epochs vector
t0_gnss = 16100; 
epochs_full = [t0_gnss; epochs];
N_full = length(epochs_full);

% Variable inizialization: estimated vector
Xest = zeros(N_full, 7);
Xest(1,:) = X0;

% Variable inizialization: covariance matrix
Pcov = zeros(7,7,N_full);
Pcov(:,:,1)= P0;

% Variable inizialization: trace of the covariance matrix
traceP = zeros(N_full,1);
traceP(1) = trace(P0);

% Variable inizialization: vectors containing the residuals for pseudo-range
% and pseudo-range rate. The first column contains the pre-fit residuals
% and the second contains the post-fit residuals
resRange   = zeros(length(epochs)  , 2 );
resRate    = zeros(length(epochs) , 2 );

% options
Tol0=1e-13;
Tol1=1e-13;
options = odeset('RelTol',Tol0,'AbsTol',Tol1);


for i = 2:length(epochs_full)

    % Integration of the 7 state equations + 49 STM equations
    X_in = Xest(i-1, :);
    phi = reshape(eye(7),1,49);
    [t, Y] = ode113(@Model_and_transition,[epochs_full(i-1) epochs_full(i)], [X_in phi], options);    
    
    % Extraction of the state and the STM at t(i)
    X_pred = Y(end, 1:7)';                   
    Phi = reshape(Y(end, 8:56), 7, 7);    % STM from t_{i-1} to t_i

    % Propagation of the covariance matrix
    Pcov(:,:,i) = Phi*Pcov(:,:,i-1)*Phi';

    % Compute the Computed Observables
    cmpObs = CmpObsSet(X_pred',epochs_full(i));
    cmp_range = cmpObs(1);
    cmp_rate  = cmpObs(2);

    % Compute pre-fit residuals
    res_range = range_obs(i-1) - cmp_range;
    res_rate = rate_obs(i-1) - cmp_rate;

    % Define the observables deviation vector
    y = [res_range; res_rate];

    % Store the residuals pre-fit
    resRange(i-1,1) = res_range;
    resRate(i-1,1) = res_rate;

    % Compute the H_tilde matrix
    [xg, yg, vxg, vyg] = gnss(propagate(epochs_full(i)));
    Ht = Htilde(X_pred, [xg, yg, vxg, vyg] );

    % Kalman gain matrix
    K = Pcov(:,:,i)*Ht.'*inv( Ht * Pcov(:,:,i) * Ht.' + R);

    % Correction
    xhat = K * y;
    Pcov(:,:,i) = (eye(7) - K * Ht ) * Pcov(:,:,i);
    traceP(i) = trace(Pcov(:,:,i));

    % Correct the state
    Xest(i,:) = (X_pred+xhat).';

    % Compute the Computed Observables with the corrected state
    cmpObs = CmpObsSet(Xest(i,:),epochs_full(i));
    cmp_range = cmpObs(1);
    cmp_rate  = cmpObs(2);

    % Compute post-fit residuals
    res_range = range_obs(i-1) - cmp_range;
    res_rate = rate_obs(i-1) - cmp_rate;

    % Store the residuals pre-fit
    resRange(i-1,2) = res_range;
    resRate(i-1,2) = res_rate;


end

%% PLOTS C1: ESTIMATED STATE VS TIME AND RESIDUALS VS TIME

% Estimated state components vs time (row 1 of Xest = a priori at t0_gnss)
state_labels = {'x [km]', 'y [km]', 'v_x [km/s]', 'v_y [km/s]', ...
                '\rho_0 [kg/km^3]', '\delta t [s]', 'd(\delta t)/dt [s/s]'};
state_titles = {'x', 'y', 'v_x', 'v_y', '\rho_0', 'Clock bias', 'Clock drift'};

figure(2)
for k = 1:7
    subplot(4,2,k)
    plot(epochs_full, Xest(:,k), 'b-')
    grid on
    title(state_titles{k})
    xlabel('Time [s]')
    ylabel(state_labels{k})
end
sgtitle('EKF estimated state components')

% Residuals vs time (converted to m and m/s for readability)
figure(3)
subplot(2,1,1)
plot(epochs, resRange(:,1)*1e3, 'r.-'); hold on
plot(epochs, resRange(:,2)*1e3, 'k.-')
yline( 0.1, 'g--'); yline(-0.1, 'g--')      % +/- 10 cm (1-sigma of the measurement)
grid on
title('Pseudo-range residuals')
xlabel('Time [s]')
ylabel('Range residual [m]')
ylim([-1 1])
legend('Pre-fit', 'Post-fit', '\pm 1\sigma (10 cm)', 'Location', 'best')

subplot(2,1,2)
plot(epochs, resRate(:,1)*1e3, 'r.-'); hold on
plot(epochs, resRate(:,2)*1e3, 'k.-')
yline( 0.01, 'g--'); yline(-0.01, 'g--')    % +/- 1 cm/s (1-sigma of the measurement)
grid on
title('Pseudo-range rate residuals')
xlabel('Time [s]')
ylabel('Range-rate residual [m/s]')
ylim([-0.1 0.1])
legend('Pre-fit', 'Post-fit', '\pm 1\sigma (1 cm/s)', 'Location', 'best')

%% C2: LINEAR REGRESSION OF THE CLOCK BIAS AND ESTIMATE OF dt0 AT t0_GNSS

% Clock bias timeline estimated by the EKF (only filtered epochs: row 1 is
% the a priori at t0_gnss, not an estimate from the data)
t_fit  = epochs;               % [s]
dt_fit = Xest(2:end, 6);       % [s]

% Linear model: dt(t) = dt0 + dt_dot*(t - t0_gnss)
% (time referred to t0_gnss, so the intercept is directly dt0)
A_fit = [ones(N_points,1), (t_fit - t0_gnss)];
coef  = A_fit \ dt_fit;
dt0_reg    = coef(1);          % bias at t0_gnss [s]
dtdot_reg  = coef(2);          % drift [s/s]

% Uncertainty of the regression parameters
res_fit  = dt_fit - A_fit*coef;
s2       = (res_fit.'*res_fit)/(N_points - 2);
cov_fit  = s2*inv(A_fit.'*A_fit);
sig_dt0_reg   = sqrt(cov_fit(1,1));
sig_dtdot_reg = sqrt(cov_fit(2,2));

% Print results
fprintf('\n===========  LINEAR REGRESSION OF THE CLOCK BIAS ===========\n');
fprintf('dt0 at t0_gnss = %.0f s : %.6e s  (+/- %.2e s)  ->  c*dt0 = %.4f km\n', ...
    t0_gnss, dt0_reg, sig_dt0_reg, c*dt0_reg);
fprintf('Drift (slope)           : %.6e s/s (+/- %.2e s/s)\n', dtdot_reg, sig_dtdot_reg);
fprintf('EKF final drift estimate: %.6e s/s\n', Xest(end,7));
fprintf('A priori dt0 / drift    : %.6e s / %.6e s/s\n', dt0, dt_dot0);
fprintf('Regression residual RMS : %.3e s\n', sqrt(mean(res_fit.^2)));
fprintf('================================================================\n');

% Plot
t_line = linspace(t0_gnss, epochs(end), 200);
dt_line = dt0_reg + dtdot_reg*(t_line - t0_gnss);

figure(4)
plot(t_fit, dt_fit, 'b.', 'MarkerSize', 8); hold on
plot(t_line, dt_line, 'r-')
plot(t0_gnss, dt0_reg, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 9)
grid on
box on
title('Clock bias: EKF estimate and linear regression')
xlabel('Time [s]')
ylabel('Clock bias \delta t [s]')
legend('EKF estimate', 'Linear regression', ...
    sprintf('\\delta t_0 = %.3e s', dt0_reg), 'Location', 'northwest')


%% PLOTS C3: TRACE OF THE COVARIANCE MATRIX VS TIME

% Prediction-only covariance between t0_gnss and the first measurement
t_pred = linspace(t0_gnss, epochs(1), 200)';
[tt, Yp] = ode113(@Model_and_transition, t_pred, [X0.' reshape(eye(7),1,49)], options);

traceP_pred = zeros(length(tt),1);
for j = 1:length(tt)
    Phi_j = reshape(Yp(j,8:56), 7, 7);          % Phi(t_j, t0_gnss)
    traceP_pred(j) = trace(Phi_j*P0*Phi_j.');
end

figure(5)
semilogy(tt, traceP_pred, 'Color', [0.2 0.4 0.9], 'LineWidth', 2); hold on
semilogy(epochs_full(2:end), traceP(2:end), 'Color', [0.85 0.1 0.1], 'LineWidth', 2)
xline(epochs(1), '--k')
grid on; box on
title('Trace of the covariance matrix')
xlabel('Time [s]')
ylabel('tr(P) [mixed units]')
legend('Prediction only (no measurements)', 'EKF (after updates)', 'Location', 'best')
