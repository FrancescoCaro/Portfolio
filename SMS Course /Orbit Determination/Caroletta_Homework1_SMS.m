clear all
close all
clc
format long e
set(groot,'defaultLineLineWidth',2)
set(groot,'DefaultAxesFontSize',15)
observables = 'both';     % This entry can be changed to 'range' or 'rate' if we want to use only the range or the range rate observables for the process

% PARAMETERS OF THE SYSTEM
global AU c Phi_1AU theta A_eff m_MPO d_mercury_sun d_earth_sun R_earth R_mercury dx dy; %#ok<GVMIS>
% Astronomical and Physical Constants
AU = 149597870.70;      % Astronomical Unit [km]
c  = 299792.458;        % Speed of light [km/s]

% Mission and Spacecraft Parameters
Phi_1AU = 1371;         % Solar flux at 1 AU [W/m^2]
A_eff   = 5.2e-6;          % Spacecraft total effective area [km^2]
m_MPO   = 1107.25;      % MPO mass [kg]

% Planetary Parameters
% Distances from the Sun
d_mercury_sun = 0.41*AU;   % Distance between Mercury and the Sun [km]
d_earth_sun   = 1.0*AU;    % Distance between the Earth and the Sun [km]

% Radii
R_earth   = 6378.14;    % Earth radius [km]
R_mercury = 2440.0;     % Mercury radius [km]

theta = 18;   % Angle between the heliocentric position of the Earth and the heliocentric position of Mercury [deg]

dx = d_mercury_sun * cosd(theta);  
dy = d_mercury_sun * sind(theta);



%% OBSERVABLES LOADING
obs_file = load('observables.txt');
epochs = obs_file(:,1);     % Time tags: first column of observable file
range_obs = obs_file(:,2);  % Range observables: second column of observable file [km]
rate_obs = obs_file(:,3);   % Range rate observables: third column of observabe file [km/s]
N_points = length(epochs);

% Observables plot
figure(1)
subplot(2,1,1)
plot(epochs, range_obs, 'k*-')
title('Range observed observables')
xlabel('Time [s]')
ylabel('Range [km]')

subplot(2,1,2)
plot(epochs, rate_obs, 'r*-')
title('Range-rate observed observables')
xlabel('Time [s]')
ylabel('Range rate [km/s]')

%% DEFINITION OF THE FUNCTIONS OF THE DYNAMICAL SYSTEM

function A = dynamical_matrix(X) 
    % Compute the dynamical matrix A starting from the state vector X
    %
    % Input:
    %   - X: state vector
    %
    % Output:
    %   - A: dynamical matrix

    global Phi_1AU c theta AU m_MPO A_eff dx dy 

    
    r = sqrt(X(1)^2 + X(2)^2);        % Distance from the spacecraft to Mercury's center [km]
    Rmpo = sqrt((X(1)+dx)^2+(X(2)+dy)^2);  % Distance from the spacecraft to the Sun [km]
    K = (Phi_1AU/c) * (AU^2) * (A_eff/m_MPO) * (1 + X(6));        % Acceleration magnitude per unit of (1 + Cs)
    common_den = -2 * K / (Rmpo^4);      % Common denominator for SRP partial derivatives with respect to position

    % intialization
    A = zeros(6,6);

    % first row
    A(1,3) = 1;

    % second row
    A(2,4) = 1;

    % third row
    A(3,1) = - X(5)/r^3 + 3*X(5)*X(1)^2/r^5 + common_den * (X(1) + dx) * cosd(theta);
    A(3,2)= 3*X(5)*X(1)*X(2)/r^5 + common_den * (X(2) + dy) * cosd(theta);
    A(3,5) = - X(1)/r^3;
    A(3,6) = (Phi_1AU/c) *(AU/Rmpo)^2 * (A_eff/m_MPO) * cosd(theta) ;

    % fourht row
    A(4,1) = 3*X(5)*X(1)*X(2)/r^5 + common_den * (X(1) + dx) * sind(theta);
    A(4,2)= - X(5)/r^3 + 3*X(5)*X(2)^2/r^5 + common_den * (X(2) + dy) * sind(theta) ;
    A(4,5) = - X(2)/r^3;
    A(4,6) = (Phi_1AU/c) *(AU/Rmpo)^2 * (A_eff/m_MPO) * sind(theta) ;
end



function dX = Model_and_transition(~,X) 
    % ODE function to integrate both trajectory and state-transition-matrix
    global AU c Phi_1AU theta A_eff m_MPO dx dy;

    % initialization
    dX = zeros(42,1); % 6->state + 36->state transition matrix element

    % state integration
    dX(1) = X(3);
    dX(2) = X(4);
    dX(3) = -X(5)*X(1)/(sqrt(X(1)^2+X(2)^2))^3 + Phi_1AU/c *(AU/sqrt((X(1)+dx)^2+(X(2)+dy)^2))^2 * A_eff/m_MPO * (1+X(6)) * cosd(theta) ;
    dX(4)= -X(5)*X(2)/(sqrt(X(1)^2+X(2)^2))^3 + Phi_1AU/c *(AU/sqrt((X(1)+dx)^2+(X(2)+dy)^2))^2 * A_eff/m_MPO * (1+X(6)) * sind(theta);
    
    % STM integration
    phi = X(7:42);
    PHI = reshape(phi,6,6);
    A   = dynamical_matrix(X);
    dphi = A*PHI;
    dphi = reshape(dphi,1,36);
    
    dX(7:42) = dphi;

end

function Obs = CmpObsSet(X) 
    % Compute the range and range rate observable.
    %
    % Input:
    %   - X: [Nx6] array containing the state vector at the different
    %              epochs
    %
    % Output:
    %   - Obs: array [Nx2] containing range and range-rate

    global AU  theta  R_earth  ; %#ok<GVMIS>

    r_merc_sun = [0.41*AU*cosd(theta),0.41*AU*sind(theta)];
    r_station_earth = [AU + R_earth * cosd(120), R_earth * sind(120) ];
    r_merc_stat =   r_merc_sun - r_station_earth ;     

    % initialization
    Obs = zeros(size(X,1),2);
    x = X(:,1);
    y = X(:,2);
    vx = X(:,3);
    vy = X(:,4);

    % computation
    for i=1:size(X,1)
        range = sqrt((x(i)+r_merc_stat(1))^2+(y(i)+r_merc_stat(2))^2);
        rate = ((x(i)+r_merc_stat(1))*vx(i) + (y(i)+r_merc_stat(2))*vy(i)) / range;
        Obs(i,1) = range;
        Obs(i,2) = rate;
    end

end


function H = Htilde(Xt, Obst)        
    % Compute the range and range-rate mapping matrix.
    %
    % Input:
    %   Xt   - [1x6] or [6x1] state at the reference epoch [x, y, vx, vy, GM, Cs]
    %   Obst - [1x2] or [2x1] observable at the reference epoch [range, rate]
    %
    % Output:
    %   H    - [2x6] mapping matrix
    global AU   theta     R_earth   ; %#ok<GVMIS>
    % Initialization
    H = zeros(2,6);

    Dx = 0.41*AU*cosd(theta) -  (AU+R_earth*cosd(120));
    Dy = 0.41*AU*sind(theta)- R_earth*sind(120);

    x = Xt(1);
    y = Xt(2);
    vx = Xt(3);
    vy = Xt(4);
    range_x = x+Dx;
    range_y = y+Dy;
    range = Obst(1);
    rate = Obst(2);

    % Computation
    % First row
    H(1,1) = range_x/range;
    H(1,2) = range_y/range;
    % Second row
    H(2,1) = (vx-rate*(range_x/range))/range;
    H(2,2) = (vy-rate*(range_y/range))/range;
    H(2,3) = range_x/range;
    H(2,4) = range_y/range;
end



%% FILTER CODE
% Filter parameters
iterations = 5;  % Number of the filter iterations
X0 = [-2840.98,761.22,-0.853853,-2.802925,22032.55,0.083];
dX0_apriori = zeros(6,1);
sig_range = 1e-3;%1.046e-4;
sig_rate  = 1e-7% 1.666e-8;
range_bias = 0.0;
epochs_full = [0;epochs];

% weight definition
w_range = (1/sig_range)^2;
w_rate  = (1/sig_rate)^2;

% a priori covariance matrix
sig0 = [0.1, 0.1, 1e-5, 1e-5, 0.5, 0.01];
P0_inv = diag(1 ./ sig0.^2);

% initialization
resRange   = zeros(length(epochs)  , iterations);
resRate    = zeros(length(epochs)  , iterations);
residual   = zeros(length(epochs)*2, iterations);
iterDelta  = zeros(6, iterations);
iterState  = zeros(6, iterations);
cov        = zeros(6, 6, iterations);


% intial condition
iterState(:,1) = X0;


% options
Tol0=1e-13;
Tol1=1e-13;
options = odeset('RelTol',Tol0,'AbsTol',Tol1);

for i = 1:iterations

    % trajectory and state transition matrix computation
    X0 = iterState(:,i);
    PHI=eye(6);
    phi=reshape(PHI,1,36);

    % integrate to obtain SMS(t) and X*(t)
    [t,X]=ode113(@Model_and_transition,epochs_full,[X0',phi],options);

    % computed observables
    cmpObs    = CmpObsSet(X(2:end,1:6));
    cmp_range = cmpObs(:,1);
    cmp_rate  = cmpObs(:,2);

    % pre-fit residuals
    res_range = range_obs - cmp_range;    
    res_rate  = rate_obs  - cmp_rate;
 
    % store residual history
    resRange(:,i) = res_range;
    resRate(:,i)  = res_rate;

    % observation deviation vector
    switch observables
        case 'both'
            y = [res_range; res_rate];
        case 'range'
            y = res_range;
        case 'rate'
            y = res_rate;
        otherwise
            error('observables must be: both, range, or rate');
    end
    residual(1:length(y),i) = y;

    % generation of the H matrix
    H = zeros(length(y), 6);
    N_range = length(res_range);

    for j = 1:length(epochs)
        % state transition matrix phi(tj,t0)
        PHI = X(j+1,7:42);
        phi = reshape(PHI,6,6);

        % H-tilde at epoch tj
        Htj = Htilde(X(j+1,:), cmpObs(j,:));

        % mapping to initial epoch
        Ht0 = Htj * phi;

        % fill final mapping matrix
        switch observables
            case 'both'
                H(j,:)         = Ht0(1,:);
                H(N_range+j,:) = Ht0(2,:);
            case 'range'
                H(j,:) = Ht0(1,:);
            case 'rate'
                H(j,:) = Ht0(2,:);
            otherwise
                error('observables must be: both, range, or rate');
        end
    end

    % weight matrix
    W = zeros(length(y), length(y));
    switch observables
        case 'both'
            W(1:N_range, 1:N_range) = eye(N_range) * w_range;
            W(N_range+1:end, N_range+1:end) = eye(N_range) * w_rate;
        case 'range'
            W = eye(length(y)) * w_range;
        case 'rate'
            W = eye(length(y)) * w_rate;
        otherwise
            error('observables must be: both, range, or rate');
    end

    % filter update
    x_hat = (H' * W * H + P0_inv) \ (H' * W * y + P0_inv * dX0_apriori);
    cov(:,:,i) = (H' * W * H + P0_inv)\ eye(6);
    iterDelta(:,i) = x_hat;

    if i < iterations
        iterState(:,i+1) = iterState(:,i) + x_hat;
        dX0_apriori = dX0_apriori - x_hat;
    end
end

%% CONVERGENCE CHECK

% print of the table at each iteration
for i=1:iterations
    fprintf('\n*** Iteration: %d ***\n',i)
    fprintf('%-12s | %-12s | %-12s | %-12s | %-12s |\n','PARAMETER', 'NOMINAL', 'CORRECTED', 'DELTA', 'COMPUTED');

    % extract the values: x0
    dx0 = iterDelta(1,i);
    x0  = iterState(1,i);
    x0N = iterState(1,1);
    x0_sig = sqrt(cov(1,1,i));
    IterDeltaV0 = iterDelta(2,i);
    fprintf('%-12s | %12.5e | %12.5e | %12.5e | %12.5e | \n' ,'x0', x0N, x0, dx0, x0_sig)

    % extract the values: y0
    dy0 = iterDelta(2,i);
    y0  = iterState(2,i);
    y0N = iterState(2,1);
    y0_sig = sqrt(cov(2,2,i));
    fprintf('%-12s | %12.5e | %12.5e | %12.5e | %12.5e | \n' ,'y0', y0N, y0, dy0, y0_sig)

     % extract the values: vx0
    dvx0 = iterDelta(3,i);
    vx0  = iterState(3,i);
    vx0N = iterState(3,1);
    vx0_sig = sqrt(cov(3,3,i));
    fprintf('%-12s | %12.5e | %12.5e | %12.5e | %12.5e | \n' ,'vx0', vx0N, vx0, dvx0, vx0_sig)

     % extract the values: vy0
    dvy0 = iterDelta(4,i);
    vy0  = iterState(4,i);
    vy0N = iterState(4,1);
    vy0_sig = sqrt(cov(4,4,i));
    fprintf('%-12s | %12.5e | %12.5e | %12.5e | %12.5e | \n' ,'vy0', vy0N, vy0, dvy0, vy0_sig)

     % extract the values: GM0
    dGM0 = iterDelta(5,i);
    GM0  = iterState(5,i);
    GM0N = iterState(5,1);
    GM0_sig = sqrt(cov(5,5,i));
    fprintf('%-12s | %12.5e | %12.5e | %12.5e | %12.5e | \n' ,'GM0', GM0N, GM0, dGM0, GM0_sig)

     % extract the values: Cs0
    dCs0 = iterDelta(6,i);
    Cs0  = iterState(6,i);
    Cs0N = iterState(6,1);
    Cs0_sig = sqrt(cov(6,6,i));
    fprintf('%-12s | %12.5e | %12.5e | %12.5e | %12.5e | \n' ,'Cs0', Cs0N, Cs0, dCs0, Cs0_sig)
end


% take the state at the different iterations
iterVec=1:1:iterations;
iter_x = iterState(1,:);
iter_y = iterState(2,:);
iter_vx   = iterState(3,:);
iter_vy   = iterState(4,:);
iter_GM   = iterState(5,:);
iter_Cs   = iterState(6,:);

% 3-sigma +/- for x
plus3sigma_x = iter_x + 3*sqrt(squeeze(cov(1,1,:)))';
minus3sigma_x = iter_x - 3*sqrt(squeeze(cov(1,1,:)))';

% 3-sigma +/- for y
plus3sigma_y = iter_y + 3*sqrt(squeeze(cov(2,2,:)))';
minus3sigma_y = iter_y - 3*sqrt(squeeze(cov(2,2,:)))';

% 3-sigma +/- for vx
plus3sigma_vx = iter_vx + 3*sqrt(squeeze(cov(3,3,:)))';
minus3sigma_vx = iter_vx - 3*sqrt(squeeze(cov(3,3,:)))';

% 3-sigma +/- for vy
plus3sigma_vy = iter_vy + 3*sqrt(squeeze(cov(4,4,:)))';
minus3sigma_vy = iter_vy - 3*sqrt(squeeze(cov(4,4,:)))';

% 3-sigma +/- for GM
plus3sigma_GM = iter_GM + 3*sqrt(squeeze(cov(5,5,:)))';
minus3sigma_GM = iter_GM - 3*sqrt(squeeze(cov(5,5,:)))';

% 3-sigma +/- for Cs
plus3sigma_Cs = iter_Cs + 3*sqrt(squeeze(cov(6,6,:)))';
minus3sigma_Cs = iter_Cs - 3*sqrt(squeeze(cov(6,6,:)))';

figure(2)
state_names = {'x [km]', 'y [km]', 'vx [km/s]', 'vy [km/s]', 'GM [km^3/s^2]', 'Cs'};
iter_states = {iter_x, iter_y, iter_vx, iter_vy, iter_GM, iter_Cs};
plus3 = {plus3sigma_x, plus3sigma_y, plus3sigma_vx, plus3sigma_vy, plus3sigma_GM, plus3sigma_Cs};
minus3 = {minus3sigma_x, minus3sigma_y, minus3sigma_vx, minus3sigma_vy, minus3sigma_GM, minus3sigma_Cs};

for k = 1:6
    subplot(3,2,k)
    hold on
    patch([iterVec fliplr(iterVec)], [minus3{k} fliplr(plus3{k})], 'g', ...
          'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(iterVec, iter_states{k}, 'r*-', 'LineWidth', 1.5);
    title(['Convergence of ' state_names{k}]);
    xlabel('Iterations');
    ylabel(state_names{k});
    grid on;
    if k == 1
        legend('3\sigma Uncertainty', 'Estimated value', 'Location', 'best');
    end
    hold off
end
sgtitle('Convergence of State Vector Components');





%% RESIDUALS PLOT
for i=1:iterations
    res = residual(:,i);
    % extract the residuals
    res_range = res(1:N_points);
    res_rate  = res(N_points+1:end);
    % ecomute the mean
    mean_range = mean(res_range);
    mean_rate  = mean(res_rate);
    % compute RMS
    RMS_range = sqrt(res_range'*res_range/length(res_range));
    RMS_rate  = sqrt(res_rate'*res_rate/length(res_rate));
    % compute SOS
    SOS_range  = res_range'*w_range*res_range;
    SOS_rate   = res_rate'*w_rate*res_rate;
    SOS = res_range'*w_range*res_range + res_rate'*w_rate*res_rate;


    switch observables
        case 'both'
            fprintf('\n*** Iteration: %d',i);
            figure('DefaultAxesFontSize',12);
            string=sprintf('Range residuals: Iter %i\n M: %1.2e RMS: %1.3e  SOS: %1.3e',[i, mean_range, RMS_range, SOS_range]);
            plot(epochs, res_range, '+')
            title(string)
            xlabel('time [s]')
            ylabel('[km]')
            
            figure('DefaultAxesFontSize',12);
            string=sprintf('Range-rate residuals: Iter %i\n M: %1.2e RMS: %1.3e  SOS: %1.3e',[i, mean_rate, RMS_rate, SOS_rate]);
            plot(epochs, res_rate,'+')
            title(string)
            xlabel('time [s]')
            ylabel('[km/s]')  

        case 'range'
            fprintf('\n*** Iteration: %d',i);
            figure('DefaultAxesFontSize',12);
            string=sprintf('Range residuals: Iter %i\n M: %1.2e RMS: %1.3e  SOS: %1.3e',[i, mean_range, RMS_range, SOS_range]);
            plot(epochs, res_range, '+')
            title(string)
            xlabel('time [s]')
            ylabel('[km]')

        case 'rate'
            res_rate = res(1:N_points);
            % compute the mean
            mean_rate  = mean(res_rate);
            % compute SOS
            SOS_rate   = res_rate'*w_rate*res_rate;
            % compute RMS
            RMS_rate  = sqrt(res_rate'*res_rate/length(res_rate));
            fprintf('\n*** Iteration: %d',i);
            figure('DefaultAxesFontSize',12);
            string=sprintf('Range-rate residuals: Iter %i\n M: %1.2e RMS: %1.3e  SOS: %1.3e',[i, mean_rate, RMS_rate, SOS_rate]);
            plot(epochs, res_rate,'+')
            title(string)
            xlabel('time [s]')
            ylabel('[km/s]')

        otherwise
            error('observables must be: both, range, or rate');
    end
end

%% FINAL ESTIMATE PRINTOUT 
% Extract the final state vector at the last iteration (t = 0)
final_state = iterState(:, iterations);

% Extract the final covariance matrix
final_cov = cov(:,:,iterations);

% Compute the 1-sigma uncertainties (square root of the diagonal elements)
final_sigma = sqrt(diag(final_cov));

% Define parameter names and units for formatting
param_names = {'x0 [km]', 'y0 [km]', 'vx0 [km/s]', 'vy0 [km/s]', 'GM0 [km^3/s^2]','Cs0 '};

% Print the final results table
fprintf('\n======================================================\n');
fprintf('   FINAL STATE ESTIMATE AT REFERENCE EPOCH (t = 0)    \n');
fprintf('======================================================\n');
fprintf('%-20s | %-15s | %-15s\n', 'PARAMETER', 'ESTIMATED VALUE', '1-SIGMA UNCERT.');
fprintf('------------------------------------------------------\n');

% Loop through the 6 state components and print values
for k = 1:6
    fprintf('%-20s | %15.6e | %15.6e\n', param_names{k}, final_state(k), final_sigma(k));
end
fprintf('======================================================\n');

%% OBSERVABLES COMPARISON  
% Extract the final covariance matrix and 1-sigma uncertainties
final_cov = cov(:,:,iterations);
final_sigma = sqrt(diag(final_cov));

fprintf('\n======================================================\n');
fprintf('   UNCERTAINTY ANALYSIS FOR OBSERVABLE SELECTION      \n');
fprintf('======================================================\n');

switch observables
    case 'range'
        fprintf('Current observable: RANGE ONLY\n');
        fprintf('%-12s | %-20s\n', 'PARAMETER', 'UNCERTAINTY');
        fprintf('------------------------------------------------------\n');
        fprintf('%-12s | %15.6e [km]\n'      , 'x0', final_sigma(1));
        fprintf('%-12s | %15.6e [km]\n'      , 'y0', final_sigma(2));
        fprintf('%-12s | %15.6e [km/s]\n'    , 'vx0', final_sigma(3));
        fprintf('%-12s | %15.6e [km/s]\n'    , 'vy0', final_sigma(4));
        fprintf('%-12s | %15.6e [km^3/s^2]\n' , 'GM0', final_sigma(5));
        fprintf('%-12s | %15.6e \n' , 'Cs0', final_sigma(6));
        
    case 'rate'
        fprintf('Current observable: RANGE-RATE ONLY\n');
        fprintf('%-12s | %-20s\n', 'PARAMETER', 'UNCERTAINTY');
        fprintf('------------------------------------------------------\n');
        fprintf('%-12s | %15.6e [km]\n'      , 'x0', final_sigma(1));
        fprintf('%-12s | %15.6e [km]\n'      , 'y0', final_sigma(2));
        fprintf('%-12s | %15.6e [km/s]\n'    , 'vx0', final_sigma(3));
        fprintf('%-12s | %15.6e [km/s]\n'    , 'vy0', final_sigma(4));
        fprintf('%-12s | %15.6e [km^3/s^2]\n' , 'GM0', final_sigma(5));
        fprintf('%-12s | %15.6e \n' , 'Cs0', final_sigma(6));
        
    case 'both'
        fprintf('Current observable: BOTH\n');
        fprintf('NOTE: To answer Point 2 of the Homework, change the\n');
        fprintf('"observables" variable at the top of the script to\n');
        fprintf('''range'' or ''rate'' and run the code again.\n');
end
fprintf('======================================================\n');
%% ALTITUDE AND UNCERTAINTY OVER TIME

% integrate final trajectory
X0_final = iterState(:,end);
PHI = eye(6);
phi = reshape(PHI,1,36);
[~, X_final] = ode113(@Model_and_transition, epochs_full, [X0_final', phi], options);

% altitude
x_traj = X_final(2:end, 1);
y_traj = X_final(2:end, 2);
r_traj = sqrt(x_traj.^2 + y_traj.^2);
alt    = r_traj - R_mercury;

% uncertainty on altitude via error propagation
P_final  = cov(:,:,end);
sig_alt  = zeros(length(epochs), 1);

for k = 1:length(epochs)
    PHI_k = reshape(X_final(k+1, 7:42), 6, 6);
    P_k   = PHI_k * P_final * PHI_k';
    x_k   = x_traj(k);
    y_k   = y_traj(k);
    r_k   = r_traj(k);
    grad  = [x_k/r_k, y_k/r_k, 0, 0, 0, 0];
    sig_alt(k) = sqrt(grad * P_k * grad');
end

figure()
subplot(2,1,1)
plot(epochs, alt, 'b-')
xlabel('Time [s]')
ylabel('Altitude [km]')
title('Altitude over time')
grid on

subplot(2,1,2)
plot(epochs, sig_alt, 'r-')
xlabel('Time [s]')
ylabel('1\sigma altitude uncertainty [km]')
title('Altitude uncertainty over time')
grid on


%%   GAP IN THE OBSERVABLES PLOT
t_plot = linspace(0, epochs(end), 2000);
[~, X_plot] = ode113(@Model_and_transition, t_plot, [X0_final', phi], options);

x_smooth = X_plot(:, 1);
y_smooth = X_plot(:, 2);

idx_gap_plot = t_plot >= 7260 & t_plot <= 9600;

%  Plot
figure()
hold on

% Mercury
theta_circle = linspace(0, 2*pi, 100);
fill(R_mercury*cos(theta_circle), R_mercury*sin(theta_circle), ...
     [0.5 0.5 0.5], 'EdgeColor', 'none', 'DisplayName', 'Mercury');

% full trajectory 
plot(x_smooth, y_smooth, 'b-', 'DisplayName', 'Trajectory');

% highlight gap 
plot(x_smooth(idx_gap_plot), y_smooth(idx_gap_plot), ...
     'r-', 'LineWidth', 3, 'DisplayName', 'Occultation gap');

% labels
plot(0, 0, 'k+', 'MarkerSize', 15, 'DisplayName', 'Mercury CoM');

hold off
axis equal
xlabel('X [km]')
ylabel('Y [km]')
title('Spacecraft trajectory and occultation gap')
legend('Location', 'best')
grid on









