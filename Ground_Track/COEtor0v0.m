function [r0_COE,v0_COE] = COEtor0v0(a,ecc,i,raan,w,v)

%L'input "a" è in chilometri. Tutti gli altri input sono in gradi. Gli
%output sono in Km e in Km/s

mu=398600.4355067;
i=deg2rad(i);
raan=deg2rad(raan);
w=deg2rad(w);
v=deg2rad(v);

R3_RAAN = [cos(raan), sin(-raan), 0; sin(raan), cos(raan), 0; 0, 0, 1];
R1_inc = [1, 0, 0; 0, cos(i), -sin(i); 0, sin(i), cos(i)];
R3_omega = [cos(w), sin(-w), 0; sin(w), cos(w), 0; 0, 0, 1];

R = R3_RAAN*R1_inc*R3_omega;

r = a*(1-ecc^2)/(1+ecc*cos(v));
h = sqrt(a*mu*(1-ecc^2));

r0_COE = R*r*[cos(v); sin(v); 0];
v0_COE = R*[-mu*sin(v)/h; mu*(ecc+cos(v))/h; 0];





