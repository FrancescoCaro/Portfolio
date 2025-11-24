% Propagation interval (seconds since t0)
step = 30 ; % time intervall between two points on the ground track
period = 24; % orbital period in Hours [h]
tspan = 0:step:period*60*60; % the last value corresponds to the last istant seen to draw the ground track
X0 = [r0_COE' v0_COE']; % Initial state vector 

% Propagate orbit using classical two-body problem
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[t,X] = ode45('dyn_mod',tspan,X0,options);
rVec = X(:,1:3); % ECI J2000 position vector
vVec = X(:,4:6); % ECI J2000 velocity vector


r3=rVec(:,3); 
r2=rVec(:,2); 
r1=rVec(:,1); 
norm_r=[];
wearth=2*pi/86164;
for i=1:length(r3)
    norm_r(i)=norm(rVec(i,:)); % generate a matrix whose i-th component is the norm of the position vector at the i-th time step
end
norm_r=norm_r(:); % Reshape/convert the matrix norm_r into a column vector


lat=asin(r3./norm_r); % compute the latitude

lon_eci=atan2(r2./(norm_r.*cos(lat)),r1./(norm_r.*cos(lat))); % longitude in the ECI reference frame

rot_tot=wearth*(tspan)'; % Earth angular rate
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
 
nexttile;
imagesc([-180 180],[-90 90],flipud(img))
set(gca,'YDir','normal')
hold on
plot(long,lat,'o', ...
    'MarkerEdgeColor', 'r', ...  
    'MarkerFaceColor', 'r', ...   
    'MarkerSize',2,...
    'LineStyle', 'none')
hold on
plot(long(1),lat(1),'o', ... 
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

