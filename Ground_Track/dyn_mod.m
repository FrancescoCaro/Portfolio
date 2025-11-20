function dX = dyn_mod(t,X,varargin)
mu = 398600.435507; % [km^3/s^2]

rVec = X(1:3); % position vector
vVec = X(4:6); % velocity vector

% Derivative of the state vector [rVec;vVec]
r = norm(rVec);
aVec = -mu/r^3 * rVec;
dX = [vVec;aVec];
end
