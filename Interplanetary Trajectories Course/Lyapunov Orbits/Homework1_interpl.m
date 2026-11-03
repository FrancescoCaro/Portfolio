clc, clear all, close all 
format long e
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesFontSize',15)
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

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
it1 = 0;
while true
    it1 = it1+1;
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
it2 = 0;
while true
    it2 = it2+1;
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


%% NATURAL CONTINUATION

dx0   = 1e-5;     % initial step (first orbit has a very small amplitude)
dxmax = 2e-3;     % maximum step (keeps adjacent orbits close for the bisection)
dxmin = 1e-9;     % minimum step before giving up
Nmax  = 3000;     % maximum number of orbits per family
Cmin  = 3.130;    % stop when C goes below this value

% Each row of fam is [x0 v0 T C]
fam1 = natural_continuation(x01, v01, T01, xL1, dx0, dxmax, dxmin, Nmax, Cmin, options);
fam2 = natural_continuation(x02, v02, T02, xL2, dx0, dxmax, dxmin, Nmax, Cmin, options);

fprintf('L1: %d orbits, C in [%.5f, %.5f]\n', size(fam1,1), min(fam1(:,4)), max(fam1(:,4)));
fprintf('L2: %d orbits, C in [%.5f, %.5f]\n', size(fam2,1), min(fam2(:,4)), max(fam2(:,4)));

%% PAIRS WITH THE SAME JACOBI CONSTANT (bisection on x0)

Cstar = [3.1380 3.1420 3.1470];     % target values inside [3.1370, 3.1493]
nP = numel(Cstar);
orb1 = zeros(nP, 4);                % [x0 v0 T C] for L1
orb2 = zeros(nP, 4);                % [x0 v0 T C] for L2

for i = 1:nP
    orb1(i,:) = bisection_C(fam1, Cstar(i), options);
    orb2(i,:) = bisection_C(fam2, Cstar(i), options);
end

% Summary table
fprintf('\n  C*            | L1: x0, v0, T, C | L2: x0, v0, T, C\n');
for i = 1:nP
    fprintf('%.4f | %.10f %.10f %.8f %.12f | %.10f %.10f %.8f %.12f\n', ...
        Cstar(i), orb1(i,:), orb2(i,:));
end

%% VERIFICATION AND PLOTS

colors = lines(nP);

% Closure and Jacobi checks for the paired orbits
for i = 1:nP
    for j = 1:2
        if j == 1, o = orb1(i,:); else, o = orb2(i,:); end
        [closure, dC] = check_orbit(o(1:3), options);
        fprintf('C*=%.4f L%d: closure = %.2e, max|C-C0| = %.2e\n', ...
            Cstar(i), j, closure, dC);
    end
end

% Window around the Moon, covering L1, L2 and the paired orbits
xa = xL1 - 0.04;
xb = xL2 + 0.04;
ya = 0.15;

% Grid for the zero-velocity curves (ZVC): 2U(x,y) = C
[Xg, Yg] = meshgrid(linspace(xa, xb, 900), linspace(-ya, ya, 900));
Ug = U_grid(Xg, Yg);

figure(10), clf, hold on, axis equal, grid on, box on
xlabel('x [-]'), ylabel('y [-]')


hOrb = gobjects(nP,1);
for i = 1:nP
    % ZVC (dashed)
    [~, hc] = contour(Xg, Yg, 2*Ug, [Cstar(i) Cstar(i)]);
    hc.LineColor = colors(i,:);
    hc.LineStyle = '--';
    hc.LineWidth = 1.2;

    % Paired orbits: same color = same C (solid)
    hOrb(i) = plot_orbit(orb1(i,1:3), options, colors(i,:), 2);   % L1
    plot_orbit(orb2(i,1:3), options, colors(i,:), 2);             % L2
end

% Moon and libration points
plot(1-mu, 0, 'ko', 'MarkerFaceColor', [.5 .5 .5], 'MarkerSize', 9)
plot(xL1, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2)
plot(xL2, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2)
text(xL1, -0.015, 'L_1', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
text(xL2, -0.015, 'L_2', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
text(1-mu, -0.015, 'Moon', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')

% Legend with the C values (proxy handles for the line styles)
hO = plot(NaN, NaN, 'k-',  'LineWidth', 2);
hZ = plot(NaN, NaN, 'k--', 'LineWidth', 1.2);
legend([hOrb; hO; hZ], [compose('C = %.4f', Cstar(:)); {'Lyapunov orbits'}; {'ZVC'}], ...
       'Location', 'northeastoutside')

xlim([xa xb]), ylim([-ya ya])
drawnow


%% OUTPUT FOR THE REPORT
TU = sqrt(R^3/(GM_Earth+GM_Moon))/86400;    % time unit in days
fprintf('\nmu = %.15e\n', mu);
fprintf('xL1 = %.12f   xL2 = %.12f\n', xL1, xL2);
fprintf('C_L1 = %.8f   C_L2 = %.8f\n', 2*potential(xL1,0), 2*potential(xL2,0));
disp('lambda L1:'), disp(lambda1.')
disp('lambda L2:'), disp(lambda2.')
fprintf('omega1 = %.10f  T_lin1 = %.10f\n', omega1, 2*pi/omega1);
fprintf('omega2 = %.10f  T_lin2 = %.10f\n', omega2, 2*pi/omega2);

% Initial Converged Orbits
fprintf('\nL1 initial: x0=%.12f v0=%.12e T=%.12f C=%.10f (it=%d)\n', ...
    X01(1), X01(2), X01(3), jacobi_const(X01(1),X01(2)), it1);
fprintf('L2 initial: x0=%.12f v0=%.12e T=%.12f C=%.10f (it=%d)\n', ...
    X02(1), X02(2), X02(3), jacobi_const(X02(1),X02(2)), it2);

% Families
for j = 1:2
    if j==1, f = fam1; else, f = fam2; end
    fprintf('L%d: N=%d, x0 in [%.6f, %.6f], C in [%.6f, %.6f], T in [%.4f, %.4f] (%.2f-%.2f days)\n', ...
        j, size(f,1), min(f(:,1)), max(f(:,1)), min(f(:,4)), max(f(:,4)), ...
        min(f(:,3)), max(f(:,3)), min(f(:,3))*TU, max(f(:,3))*TU);
end

% Pairs with bisection C
fprintf('\nPairs:\n');
for i = 1:nP
    [~, n1] = bisection_C(fam1, Cstar(i), options);
    [~, n2] = bisection_C(fam2, Cstar(i), options);
    fprintf('C*=%.4f | dC=%.2e | bisection it L1=%d L2=%d | T1=%.4f d, T2=%.4f d\n', ...
        Cstar(i), abs(orb1(i,4)-orb2(i,4)), n1, n2, orb1(i,3)*TU, orb2(i,3)*TU);
end

%% REPORT PLOTS
cm = parula(25);

figure(11), clf, hold on, axis equal, grid on, box on
i1 = round(linspace(1,size(fam1,1),25));
i2 = round(linspace(1,size(fam2,1),25));
for j = 1:25
    plot_orbit(fam1(i1(j),1:3), options, cm(j,:), 1);
    plot_orbit(fam2(i2(j),1:3), options, cm(j,:), 1);
end
plot(1-mu,0,'ko','MarkerFaceColor',[.5 .5 .5],'MarkerSize',9)
plot([xL1 xL2],[0 0],'kx','MarkerSize',10)
xlabel('x [-]'), ylabel('y [-]')
exportgraphics(gcf,'figFamilies.png','Resolution',300)

figure(12), clf, hold on, grid on, box on
plot(fam1(:,1), fam1(:,4), 'DisplayName','L_1')
plot(fam2(:,1), fam2(:,4), 'DisplayName','L_2')
xlabel('x_0 [-]'), ylabel('C [-]'), legend
exportgraphics(gcf,'figCvsX0.png','Resolution',300)

figure(13), clf, hold on, grid on, box on
plot(fam1(:,3)*TU, fam1(:,4), 'DisplayName','L_1')
plot(fam2(:,3)*TU, fam2(:,4), 'DisplayName','L_2')
xlabel('T [days]'), ylabel('C [-]'), legend
exportgraphics(gcf,'figCvsT.png','Resolution',300)

exportgraphics(figure(10),'figPairs.png','Resolution',300)



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

    % Extract state for readability
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
    
    % Initialization 
    DF = zeros(2,3);

    % Extract state for readability
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


function C = jacobi_const(x0,v0)
    U = potential(x0,0);
    C = 2*U - v0^2;
end


% Single shooting with FIXED x0: unknowns are [v0; T]
function [v0, T, ok] = shooting_fixed_x0(x0, v0, T, options)
tolF = 1e-12;       % tolerance on the constraint vector
maxIter = 20;       % maximum number of Newton iterations
ok = false;
tau = [0 0.5];

for it = 1:maxIter
    % Reject non-physical periods (T <= 0 is also the trivial solution)
    if T <= 0 || ~isfinite(T) || ~isfinite(v0)
        return
    end

    X0_full = [x0; 0; 0; v0; reshape(eye(4,4),16,1)];

    % Integration of the state and STM until tau = 0.5
    [~, X] = ode113(@(t,X) Model_and_transition(t,X,T), tau, X0_full, options);
    Xf = X(end,:);

    % Constraint vector: y(T/2) = 0, u(T/2) = 0
    F = [Xf(2); Xf(3)];

    if norm(F, inf) < tolF
        ok = true;
        return
    end

    % STM at tau = 0.5 and Jacobian DF (2x3)
    PHI_half = reshape(Xf(5:20),4,4);
    DF = nonlin_jac(PHI_half, [x0; 0; 0; v0]);

    % x0 is fixed: drop its column, solve the 2x2 system in [v0; T]
    d = -DF(:,2:3) \ F;
    v0 = v0 + d(1);
    T  = T  + d(2);
end
end

% Natural continuation on x0 with adaptive step
function fam = natural_continuation(x0, v0, T, xL, dx0, dxmax, dxmin, Nmax, Cmin, options)
% Step direction: away from the libration point
s = sign(x0 - xL);

fam = zeros(Nmax, 4);
fam(1,:) = [x0 v0 T jacobi_const(x0, v0)];
k = 1;
dx = dx0;

while k < Nmax && fam(k,4) > Cmin
    xn = fam(k,1) + s*dx;

    % Initial guess for [v0 T]
    if k >= 2
        % Secant predictor, scaled by the ratio of the steps
        hprev = abs(fam(k,1) - fam(k-1,1));
        g = fam(k,2:3) + (fam(k,2:3) - fam(k-1,2:3)) * (dx/hprev);
    else
        % First step: v0 is proportional to (x0 - xL) in the linear regime
        g = [fam(k,2)*(xn - xL)/(fam(k,1) - xL), fam(k,3)];
    end

    [vn, Tn, ok] = shooting_fixed_x0(xn, g(1), g(2), options);

    % Sanity check: the period must stay close to the previous one
    % (avoids jumping to another family or to the trivial solution)
    if ok && (Tn <= 0 || abs(Tn - fam(k,3)) > 0.1*fam(k,3))
        ok = false;
    end

    if ~ok
        % Halve the step and retry from the same orbit
        dx = dx/2;
        if dx < dxmin
            warning('NC: step below dxmin at k = %d (x0 = %.6f)', k, fam(k,1));
            break
        end
        continue
    end

    % Accept the orbit and enlarge the step (up to dxmax)
    k = k + 1;
    fam(k,:) = [xn vn Tn jacobi_const(xn, vn)];
    dx = min(1.5*dx, dxmax);
end

fam = fam(1:k,:);
end

% Bisection on x0 to find the orbit of the family with C = Cstar
function [sol, nit] = bisection_C(fam, Cstar, options)
C = fam(:,4);

% Find two consecutive orbits bracketing Cstar
k = find((C(1:end-1) - Cstar) .* (C(2:end) - Cstar) <= 0, 1);
assert(~isempty(k), 'Cstar = %.4f is outside the family range', Cstar);

a = fam(k,:);       % [x0 v0 T C]
b = fam(k+1,:);
tolC = 1e-11;

for it = 1:100
    nit = it;       

    % Midpoint on the natural parameter
    xm = 0.5*(a(1) + b(1));

    % Guess for [v0 T]: average of the two bracketing orbits
    g = 0.5*(a(2:3) + b(2:3));

    % Recompute the orbit at xm with the shooting
    [vm, Tm, ok] = shooting_fixed_x0(xm, g(1), g(2), options);
    assert(ok, 'Shooting did not converge during bisection');

    Cm = jacobi_const(xm, vm);
    sol = [xm vm Tm Cm];

    % Stop on tolerance (or if the interval is at machine precision)
    if abs(Cm - Cstar) < tolC || abs(b(1) - a(1)) < 1e-15
        return
    end

    % Keep the sub-interval that still brackets Cstar
    if (Cm - Cstar) * (a(4) - Cstar) > 0
        a = sol;
    else
        b = sol;
    end
end
end

% Integrate one full period and plot the orbit in the (x,y) plane
function h = plot_orbit(o, options, col, lw)
    X0 = [o(1); 0; 0; o(2); reshape(eye(4),16,1)];
    [~, X] = ode113(@(t,X) Model_and_transition(t,X,o(3)), linspace(0,1,400), X0, options);
    h = plot(X(:,1), X(:,2), 'Color', col, 'LineWidth', lw);
end

% Closure error and Jacobi constant drift over one full period
function [closure, dC] = check_orbit(o, options)
    X0 = [o(1); 0; 0; o(2); reshape(eye(4),16,1)];
    [~, X] = ode113(@(t,X) Model_and_transition(t,X,o(3)), linspace(0,1,400), X0, options);
    
    % Closure: state at tau = 1 must match the initial state
    closure = norm(X(end,1:4) - X(1,1:4));
    
    % Jacobi constant along the orbit: C = 2U - (u^2 + v^2)
    Cvec = zeros(size(X,1),1);
    for i = 1:size(X,1)
        Cvec(i) = 2*potential(X(i,1), X(i,2)) - (X(i,3)^2 + X(i,4)^2);
    end
    dC = max(abs(Cvec - Cvec(1)));
end

% Vectorized potential U(x,y) for matrix inputs (used for the ZVC)
function U = U_grid(x, y)
    global mu
    r1 = sqrt((x+mu).^2 + y.^2);
    r2 = sqrt((x+mu-1).^2 + y.^2);
    U = (x.^2 + y.^2)/2 + (1-mu)./r1 + mu./r2;
    U(~isfinite(U)) = NaN;
end