clear all
close all
clc
format long e
set(groot,'defaultLineLineWidth',2)
set(groot,'DefaultAxesFontSize',15)
observables = 'both';     % This entry can be changed to 'range' or 'rate' if we want to use only the range or the range rate observables for the process



%% PARAMETERS OF THE SYSTEM 
global gm Re A_eff m Cd z0 H_ref theta 

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
    % State vector: X = [x, y, vx, vy, rho0]'
    
    global gm Re A_eff m Cd z0 H_ref 
    
    % Extract state variables
    x = X(1); 
    y = X(2); 
    vx = X(3); 
    vy = X(4); 
    rho0 = X(5);
    
    % Compute derived quantities
    r = sqrt(x^2 + y^2);
    V = sqrt(vx^2 + vy^2);
    z = r - Re;
    
    % Atmospheric model
    exp_term = exp(-(z - z0)/H_ref);
    rho = rho0 * exp_term;
    k = Cd * A_eff / (2*m); 
    
    % Initialize 5x5 matrix
    A = zeros(5,5);
    
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
   
end



function dX = Model_and_transition(~,X) 
    % ODE function to integrate both trajectory and state-transition-matrix
    % State vector: X = [x, y, vx, vy, rho0]'
    
    global gm Re A_eff m Cd z0 H_ref 
    
    % Initialization (5 state elements + 25 STM elements)
    dX = zeros(30,1); 
    
    % Extract state for readability
    x = X(1); 
    y = X(2); 
    vx = X(3); 
    vy = X(4); 
    rho0 = X(5);
    
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
    
    
    % STM Integration
    % STM is a 5x5 matrix -> 25 elements -> indices 6 to 30
    phi = X(6:30);
    PHI = reshape(phi, 5, 5);
    
    A = dynamical_matrix(X);
    dphi = A * PHI;
    
    dX(6:30) = reshape(dphi, 25, 1);
end


function Obs = CmpObsSet(X) 
    % Compute the range and range rate observable.
    %
    % Input:
    %   - X: [Nx5] array containing the state vector at the different
    %              epochs
    %
    % Output:
    %   - Obs: array [Nx2] containing range and range-rate

    global Re  theta     

    % initialization
    Obs = zeros(size(X,1),2);
    x = X(:,1);
    y = X(:,2);
    vx = X(:,3);
    vy = X(:,4);

    % computation
    for i=1:size(X,1)
        range = sqrt((x(i)-Re*cosd(theta))^2+(y(i)-Re*sind(theta))^2);
        rate = ((x(i)-Re*cosd(theta))*vx(i) + (y(i)-Re*sind(theta))*vy(i)) / range;
        Obs(i,1) = range;
        Obs(i,2) = rate;
    end

end

function H = Htilde(Xt, Obst)        
    % Compute the range and range-rate mapping matrix.
    %
    % Input:
    %   Xt   - [1x5] or [5x1] state at the reference epoch [x, y, vx, vy, rho0]
    %   Obst - [1x2] or [2x1] observable at the reference epoch [range, rate]
    %
    % Output:
    %   H    - [2x5] mapping matrix
    global Re  theta
    % Initialization
    H = zeros(2,5);


    x = Xt(1);
    y = Xt(2);
    vx = Xt(3);
    vy = Xt(4);
    range_x = x-Re*cosd(theta);
    range_y = y-Re*sind(theta);
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
X0 = [4.7230591e3,4.7210591e3,-6.2288709,6.2283709,1.8080000e-2];
dX0_apriori = zeros(5,1);
sig_range = 1.043e-4;   % putative 1e-5   real 1.043e-4
sig_rate  = 1.005e-5;    % putative 1e-6   real 1.005e-5
range_bias = 0.0;
epochs_full = [0;epochs];

% weight definition
w_range = (1/sig_range)^2;
w_rate  = (1/sig_rate)^2;

% a priori covariance matrix
sig0 = [0.01, 0.01, 0.001, 0.001, 0.1];
P0_inv = diag(1 ./ sig0.^2);

% initialization
resRange   = zeros(length(epochs)  , iterations);
resRate    = zeros(length(epochs)  , iterations);
residual   = zeros(length(epochs)*2, iterations);
iterDelta  = zeros(5, iterations);
iterState  = zeros(5, iterations);
cov        = zeros(5, 5, iterations);


% intial condition
iterState(:,1) = X0;


% options
Tol0=1e-13;
Tol1=1e-13;
options = odeset('RelTol',Tol0,'AbsTol',Tol1);

for i = 1:iterations

    % trajectory and state transition matrix computation
    X0 = iterState(:,i);
    PHI=eye(5);
    phi=reshape(PHI,1,25);

    % integrate to obtain SMS(t) and X*(t)
    [t,X]=ode113(@Model_and_transition,epochs_full,[X0',phi],options);

    % computed observables
    cmpObs    = CmpObsSet(X(2:end,1:5));
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
    H = zeros(length(y), 5);
    N_range = length(res_range);

    for j = 1:length(epochs)
        % state transition matrix phi(tj,t0)
        PHI = X(j+1,6:30);
        phi = reshape(PHI,5,5);

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
    cov(:,:,i) = (H' * W * H + P0_inv)\ eye(5);
    iterDelta(:,i) = x_hat;

    if i < iterations
        iterState(:,i+1) = iterState(:,i) + x_hat;
        dX0_apriori = dX0_apriori - x_hat;
    end
end


%% CONVERGENCE CHECK

% print of the table at each iteration
for i=1:iterations
    fprintf('\n*** Iteration: %d ***\n', i);
    fprintf('%-10s | %-12s | %-12s | %-12s | %-12s\n', 'PARAMETER', 'NOMINAL', 'CORRECTED', 'DELTA', 'SIGMA');
    fprintf('-----------------------------------------------------------------------\n');
    
    % extract the values: x0
    dx0 = iterDelta(1,i);
    x0  = iterState(1,i);
    x0N = iterState(1,1);
    x0_sig = sqrt(cov(1,1,i));
    fprintf('%-10s | %12.5e | %12.5e | %12.5e | %12.5e\n' ,'x0', x0N, x0, dx0, x0_sig);
    
    % extract the values: y0
    dy0 = iterDelta(2,i);
    y0  = iterState(2,i);
    y0N = iterState(2,1);
    y0_sig = sqrt(cov(2,2,i));
    fprintf('%-10s | %12.5e | %12.5e | %12.5e | %12.5e\n' ,'y0', y0N, y0, dy0, y0_sig);
    
    % extract the values: vx0
    dvx0 = iterDelta(3,i);
    vx0  = iterState(3,i);
    vx0N = iterState(3,1);
    vx0_sig = sqrt(cov(3,3,i));
    fprintf('%-10s | %12.5e | %12.5e | %12.5e | %12.5e\n' ,'vx0', vx0N, vx0, dvx0, vx0_sig);
    
    % extract the values: vy0
    dvy0 = iterDelta(4,i);
    vy0  = iterState(4,i);
    vy0N = iterState(4,1);
    vy0_sig = sqrt(cov(4,4,i));
    fprintf('%-10s | %12.5e | %12.5e | %12.5e | %12.5e\n' ,'vy0', vy0N, vy0, dvy0, vy0_sig);
    
    % extract the values: rho0
    drho0 = iterDelta(5,i);
    rho0  = iterState(5,i);
    rho0N = iterState(5,1);
    rho0_sig = sqrt(cov(5,5,i));
    fprintf('%-10s | %12.5e | %12.5e | %12.5e | %12.5e\n' ,'rho0', rho0N, rho0, drho0, rho0_sig);
end



% take the state at the different iterations
iterVec=1:1:iterations;
iter_x    = iterState(1,:);
iter_y    = iterState(2,:);
iter_vx   = iterState(3,:);
iter_vy   = iterState(4,:);
iter_rho0 = iterState(5,:); 

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

% 3-sigma +/- for rho0
plus3sigma_rho0 = iter_rho0 + 3*sqrt(squeeze(cov(5,5,:)))';
minus3sigma_rho0 = iter_rho0 - 3*sqrt(squeeze(cov(5,5,:)))';

figure(2)
state_names = {'x [km]', 'y [km]', 'v_x [km/s]', 'v_y [km/s]', '\rho(z_0) [kg/km^3]'};
iter_states = {iter_x, iter_y, iter_vx, iter_vy, iter_rho0};
plus3 = {plus3sigma_x, plus3sigma_y, plus3sigma_vx, plus3sigma_vy, plus3sigma_rho0};
minus3 = {minus3sigma_x, minus3sigma_y, minus3sigma_vx, minus3sigma_vy, minus3sigma_rho0};

for k = 1:5
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

% residuals plot
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
            % fprintf('\n*** Iteration: %d',i);
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
            %fprintf('\n*** Iteration: %d',i);
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
            %fprintf('\n*** Iteration: %d',i);
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
param_names = {'x0 [km]', 'y0 [km]', 'vx0 [km/s]', 'vy0 [km/s]', 'rho(z0) [kg/km^3]'};

% Print the final results table
fprintf('\n======================================================\n');
fprintf('   FINAL STATE ESTIMATE AT REFERENCE EPOCH (t = 0)    \n');
fprintf('======================================================\n');
fprintf('%-20s | %-15s | %-15s\n', 'PARAMETER', 'ESTIMATED VALUE', '1-SIGMA UNCERT.');
fprintf('------------------------------------------------------\n');

% Loop through the 5 state components and print values
for k = 1:5
    fprintf('%-20s | %15.6e | %15.6e\n', param_names{k}, final_state(k), final_sigma(k));
end
fprintf('======================================================\n');



%% OBSERVABLES COMPARISON (Point 2) 
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
        fprintf('%-12s | %15.6e [kg/km^3]\n' , 'rho(z0)', final_sigma(5));
        
    case 'rate'
        fprintf('Current observable: RANGE-RATE ONLY\n');
        fprintf('%-12s | %-20s\n', 'PARAMETER', 'UNCERTAINTY');
        fprintf('------------------------------------------------------\n');
        fprintf('%-12s | %15.6e [km]\n'      , 'x0', final_sigma(1));
        fprintf('%-12s | %15.6e [km]\n'      , 'y0', final_sigma(2));
        fprintf('%-12s | %15.6e [km/s]\n'    , 'vx0', final_sigma(3));
        fprintf('%-12s | %15.6e [km/s]\n'    , 'vy0', final_sigma(4));
        fprintf('%-12s | %15.6e [kg/km^3]\n' , 'rho(z0)', final_sigma(5));
        
    case 'both'
        fprintf('Current observable: BOTH\n');
        fprintf('NOTE: To answer Point 2 of the Challenge, change the\n');
        fprintf('"observables" variable at the top of the script to\n');
        fprintf('''range'' or ''rate'' and run the code again.\n');
end
fprintf('======================================================\n');