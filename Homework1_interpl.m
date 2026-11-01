clc, clear all, close all 
format long e
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesFontSize',15)

%% PARAMETERS DEFINITION

% Primary bodies constants
GM_Earth = 3.98600436e5; % Gravitational parameter of the Earth [km^3/sec^2]
GM_Moon = 4.902799e3; % Gravitational parameter of the Moon [km^3/sec^2]
R = 384400; % Distance between the two primary bodies [km]

% CR3BP constants
global mu
mu = GM_Moon/(GM_Moon + GM_Earth); 
alpha0 = 1e-4;

% ode113 options
options = odeset('RelTol',3e-14, 'AbsTol', 1e-14);

%% LINEARIZED MOTION SOLUTION

% Lagrange points L1 and L2 computation
Ux_fun = @(x) x - (1-mu)*(x+mu)/abs(x+mu)^3 - mu*(x+mu-1)/abs(x+mu-1)^3;
int1 = [0.5 0.95];
int2 = [1.0 1.4];
xL1 = fzero(Ux_fun,int1);
xL2 = fzero(Ux_fun,int2);

% Eigenvalues and eigenvectors computation for L1
A1 = jacobian(xL1,0);
[V1, D1] = eig(A1);
lambda1 = diag(D1);
omega1 = imag(lambda1(3));
Ec1_1 = V1(:,3);

% First guess solution for L1
state0_1 = [xL1; 0; 0; 0] + alpha0*real(Ec1_1);
x01 = state0_1(1);
v01 = state0_1(4);
T01 = 2*pi/omega1;
X01 = [x01; v01; T01];

% Eigenvalues and eigenvectors computation for L2
A2 = jacobian(xL2,0);
[V2, D2] = eig(A2);
lambda2 = diag(D2);
omega2 = imag(lambda2(3));
Ec1_2 = V2(:,3);

% First guess solution for L2
state0_2 = [xL2; 0; 0; 0] + alpha0*real(Ec1_2);
x02 = state0_2(1);
v02 = state0_2(4);
T02 = 2*pi/omega2;
X02 = [x02; v02; T02];


%% FIRST ORBIT COMPUTATION

% Integration parameters
tau = [0 0.5];
eps = 1e-12;

% Orbit around L1
while true
    X01_full = [x01;0;0;v01;reshape(eye(4,4),16,1)];
    
    % Integration of the State and STM until tau = 0.5
    [t, X] = ode113(@(t,X) Model_and_transition(t,X,T01), tau, X01_full,options);
    
    % Extraction of the final state
    Xf = X(end,:);

    % Constraint check
    y_half = Xf(2);
    u_half = Xf(3);
    F = [y_half; u_half];

    if abs(F(1))<eps && (abs(F(2))<eps)
        break
    end
    
    % Extraction of the final STM and DF computation
    PHI_half = reshape(Xf(5:20),4,4);
    DF = nonlin_jac(PHI_half,[x01;0;0;v01]);
    
    % Correction of the initial state
    deltaX = - DF\F;
    X01 = X01 + deltaX;
    x01 = X01(1);
    v01 = X01(2);
    T01 = X01(3);
end



% Orbit around L2
while true
    X02_full = [x02;0;0;v02;reshape(eye(4,4),16,1)];
    
    % Integration of the State and STM until tau = 0.5
    [t, X] = ode113(@(t,X) Model_and_transition(t,X,T02), tau, X02_full,options);
    
    % Extraction of the final state
    Xf = X(end,:);

    % Constraint check
    y_half = Xf(2);
    u_half = Xf(3);
    F = [y_half; u_half];

    if abs(F(1))<eps && (abs(F(2))<eps)
        break
    end
    
    % Extraction of the final STM and DF computation
    PHI_half = reshape(Xf(5:20),4,4);
    DF = nonlin_jac(PHI_half,[x02;0;0;v02]);
    
    % Correction of the initial state
    deltaX = - DF\F;
    X02 = X02 + deltaX;
    x02 = X02(1);
    v02 = X02(2);
    T02 = X02(3);
end











%% LOCAL FUNCTIONS

% Jacobi integral function and derivatives
function [U, Ux, Uy, Uxx, Uyy, Uxy] = potential(x, y)
    global mu

    % Distances from the two primary bodies
    r1 = sqrt((x+mu)^2 + y^2);
    r2 = sqrt((x+mu-1)^2+y^2);

    % Potential 
    U = (x^2+y^2)/2 + (1-mu)/r1 + mu/r2;

    % First derivatives
    Ux = x - (1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3;
    Uy = y - (1-mu)*y/r1^3 - mu*y/r2^3;

    % Second derivatives 
    Uxx = 1 - (1 - mu)/r1^3 - mu/r2^3 ...
      + 3*(1 - mu)*(x + mu)^2/r1^5 + 3*mu*(x - 1 + mu)^2/r2^5;
    Uyy = 1 - (1 - mu)/r1^3 - mu/r2^3 ...
      + 3*(1 - mu)*y^2/r1^5 + 3*mu*y^2/r2^5;
    Uxy = 3*(1 - mu)*(x + mu)*y/r1^5 + 3*mu*(x - 1 + mu)*y/r2^5;

end

% Jacobian matrix A
function A = jacobian(x,y)
    global mu
    
    % Distances from the two primary bodies
    r1 = sqrt((x+mu)^2 + y^2);
    r2 = sqrt((x+mu-1)^2+y^2);
    
    % Partial derivatives of the potential
    Uxx = 1 - (1 - mu)/r1^3 - mu/r2^3 ...
      + 3*(1 - mu)*(x + mu)^2/r1^5 + 3*mu*(x - 1 + mu)^2/r2^5;
    Uyy = 1 - (1 - mu)/r1^3 - mu/r2^3 ...
      + 3*(1 - mu)*y^2/r1^5 + 3*mu*y^2/r2^5;
    Uxy = 3*(1 - mu)*(x + mu)*y/r1^5 + 3*mu*(x - 1 + mu)*y/r2^5;
    
    % Initialize the A matrix
    A = zeros(4,4);

    % First row
    A(1,3) = 1;

    % Second row
    A(2,4) = 1;

    % Third row
    A(3,1) = Uxx;
    A(3,2) = Uxy;
    A(3,4) = 2;

    % Fourth row
    A(4,1) = Uxy;
    A(4,2) = Uyy;
    A(4,3) = -2;

end

% State and transition matrix 
function dX = Model_and_transition(~,X,T)
    global mu
    
    % Initialization (4 state elements + 16 STM elements)
    dX = zeros(20,1);

    % Extract state for readibility
    x = X(1);
    y = X(2);
    u = X(3);
    v = X(4);

    % Distances from the two primary bodies
    r1 = sqrt((x+mu)^2 + y^2);
    r2 = sqrt((x+mu-1)^2+y^2);

    % First derivatives of the potential U
    Ux = x - (1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3;
    Uy = y - (1-mu)*y/r1^3 - mu*y/r2^3;

    % State integration
    dX(1) = T*u;
    dX(2) = T*v;
    dX(3) = (2*v + Ux)*T;
    dX(4) = (-2*u + Uy)*T;

    % STM integration
    % STM is 4x4 -> 16 elements -> indices from 5 to 20
    phi = X(5:20);
    PHI = reshape(phi, 4,4);

    A = jacobian(x,y);
    dphi =T*A*PHI;
    

    dX(5:20) = reshape(dphi, 16, 1);

end

% Jacobian matrix of the nonlinear problem
function DF = nonlin_jac(PHI, X)
    global mu
    
    % Inizialization 
    DF = zeros(2,3);

    % Extract state for readibility
    x = X(1);
    y = X(2);
    u = X(3);
    v = X(4);

    % Distances from the two primary bodies
    r1 = sqrt((x+mu)^2 + y^2);
    r2 = sqrt((x+mu-1)^2+y^2);

    % First derivatives of the potential U
    Ux = x - (1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3;
    Uy = y - (1-mu)*y/r1^3 - mu*y/r2^3;

    % State at X0
    f0 = [ u;v; (2*v + Ux); (-2*u + Uy)];


    % First row
    DF(1,1) = PHI(2,1);
    DF(1,2) = PHI(2,4);
    DF(1,3) = 0.5*PHI(2,:)*f0;

    % Second row
    DF(2,1) = PHI(3,1);
    DF(2,2) = PHI(3,4);
    DF(2,3) = 0.5*PHI(3,:)*f0;

end

