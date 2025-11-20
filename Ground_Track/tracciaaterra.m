% Propagation interval (seconds since t0)
step = 30 ; %intervallo di tempo tra un punto e un altro della traccia a terra
period = 24; % periodo orbitale in ore (h)
tspan = 0:step:period*60*60; %l'ultimo valore corrisponde al tempo passato ad osservare l'orbita (è bene di solito porlo uguale al periodo dell'orbita)
X0 = [r0_COE' v0_COE']; %VETTORI DI STATO INIZIALE (le prime tre componenti sono la posizione, le ultime tre sono la velocità)

% Propagate orbit using classical two-body problem
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[t,X] = ode45('dyn_mod',tspan,X0,options);
rVec = X(:,1:3); % ECI J2000 position vector
vVec = X(:,4:6); % ECI J2000 velocity vector


r3=rVec(:,3); %estrae la terza componente del vettore posizione
r2=rVec(:,2); %estrae la seconda componente del vettore posizione
r1=rVec(:,1); %estrae la prima componente del vettore posizione
norm_r=[];
wearth=2*pi/86164;
for i=1:length(r3)
    norm_r(i)=norm(rVec(i,:)); %genera una matrice la cui i-esima componente è la norma del vettore posizione al tempo i-esimo
end
norm_r=norm_r(:); % rende la matrice norm_r un vettore colonna


lat=asin(r3./norm_r); %calcola la latitudine

lon_eci=atan2(r2./(norm_r.*cos(lat)),r1./(norm_r.*cos(lat))); %longitudine nel sdr ECI

rot_tot=wearth*(tspan)'; %velocità di rotazione terrestre
for i=1:length(rot_tot)
    if rot_tot(i)>pi
        rot_tot(i)=rot_tot(i)-2*pi;
    elseif rot_tot(i)>2*pi
        rot_tot(i)=rot_tot(i)-2*pi;
    end
end

lon_ecef=atan2(r2./(norm_r.*cos(lat)),r1./(norm_r.*cos(lat)))-rot_tot;

for i=1:length(lon_ecef)
    if lon_ecef(i)>pi
        lon_ecef(i)=lon_ecef(i)-2*pi;
    elseif lon_ecef(i)<-pi
        lon_ecef(i)=lon_ecef(i)+2*pi;
    end
end

long=rad2deg(lon_ecef);
lat=rad2deg(lat);

file = dir('Mercator_projection_SW.*');
img = imread(file(1).name);
  % Sostituisci con il tuo file
nexttile;
imagesc([-180 180],[-90 90],flipud(img))
set(gca,'YDir','normal')
hold on
plot(long,lat,'o', ...
    'MarkerEdgeColor', 'r', ...   % bordo rosso (può anche essere 'none')
    'MarkerFaceColor', 'r', ...   % interno rosso
    'MarkerSize',2,...
    'LineStyle', 'none')
hold on
plot(long(1),lat(1),'o', ... %evidenzia in verde il punto corrispondente alla posizione iniziale
    'MarkerEdgeColor', 'g', ...   
    'MarkerFaceColor', 'g', ...  
    'MarkerSize',2,...
    'LineStyle', 'none')
xlim([-180 180]);   
ylim([-90 90]);


figure
subplot(3,1,1)
plot(tspan,lat,'o')
title('latitudine')
subplot(3,1,2)
plot(tspan,lon_eci,'o')
title('Longitudine ECI')
subplot(3,1,3)
plot(tspan,lon_ecef,'o')
title('Longitudine ECEF')
