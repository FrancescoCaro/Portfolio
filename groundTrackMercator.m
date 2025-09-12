function [lon_deg, lat_deg] = groundTrackMercator(a, e, i, Omega, omega, nu0, ...
                                                  t0, t_end, dt, mu)
% groundTrackMercator Calcola e traccia la ground track in proiezione Mercatore
%
% INPUT:
%   a, e, i, Omega, omega, nu0 : elementi orbitali classici
%   t0, t_end, dt             : tempo iniziale, finale e passo [s]
%   mu                        : parametro gravitazionale [m^3/s^2]
%
% OUTPUT:
%   lon_deg, lat_deg : vettori di longitudine e latitudine in gradi

    if nargin < 10
        mu = 3.986004418e14;  % valore di default per la Terra
    end

    % 1. Stato iniziale in ECI
    % calcola r0 e v0 da COE (funzione sotto)
    [r0, v0] = coe2rv(a, e, i, Omega, omega, nu0, mu);
    X0 = [r0; v0];

    % 2. Propagazione orbitale
    tspan = t0:dt:t_end;
    opts  = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [t, X] = ode45(@(t,X) twoBodyDyn(t,X,mu), tspan, X0, opts);

    % 3. Estrazione posizioni ECI e rotazione in ECEF
    omega_e = 2*pi/86164;  % velocità di rotazione terrestre [rad/s]
    rECI = X(:,1:3);
    lon_deg = zeros(size(t));
    lat_deg = zeros(size(t));

    for k = 1:length(t)
        % ruota di omega_e * t
        Rz = [ cos(omega_e*t(k))  sin(omega_e*t(k)) 0;
              -sin(omega_e*t(k)) cos(omega_e*t(k)) 0;
               0                  0                 1 ];
        rECEF = Rz * rECI(k,:)';
        x = rECEF(1); y = rECEF(2); z = rECEF(3);

        % 4. latitudine e longitudine
        rho = norm(rECEF);
        lat = asin(z/rho);
        lon = atan2(y,x);

        lon_deg(k) = rad2deg(lon);
        lat_deg(k) = rad2deg(lat);
    end

    % 5. Conversione in coordinate Mercatore
    Rearth = 6378137;  % raggio terrestre [m]
    xM = Rearth * deg2rad(lon_deg);
    yM = Rearth * log(tan(pi/4 + deg2rad(lat_deg)/2));

    % 6. Plot
    figure('Color','w')
    plot(xM, yM, 'b','LineWidth',1.5)
    axis equal
    grid on
    xlabel('x [m]')
    ylabel('y [m]')
    title('Ground Track in proiezione di Mercatore')
end


%% Funzione di dinamica due corpi
function dXdt = twoBodyDyn(~, X, mu)
    r = X(1:3); v = X(4:6);
    a = -mu * r / norm(r)^3;
    dXdt = [v; a];
end


%% Conversione COE → r, v in ECI
function [r, v] = coe2rv(a, e, i, Omega, omega, nu, mu)
    % calcola distanze e velocità nel PQW
    p = a*(1 - e^2);
    r_pqw = (p/(1 + e*cos(nu))) * [cos(nu); sin(nu); 0];
    v_pqw = sqrt(mu/p) * [-sin(nu); e + cos(nu); 0];

    % matrici di rotazione
    R3_Omega = [ cos(Omega), -sin(Omega), 0;
                 sin(Omega),  cos(Omega), 0;
                 0         ,  0         , 1];
    R1_i     = [1, 0         ,     0;
                0, cos(i)   , -sin(i);
                0, sin(i)   ,  cos(i)];
    R3_omega = [ cos(omega), -sin(omega), 0;
                 sin(omega),  cos(omega), 0;
                 0         ,  0         , 1];

    Q = R3_Omega * R1_i * R3_omega;
    r = Q * r_pqw;
    v = Q * v_pqw;
end
