function thrust = u2T(u, m)
% I: u - 3x1 vector of nondimensional thrust
%    m - Mass of spacecraft (kg)
% O: thrust - 3x1 vector of thrust in (kN)

D = 384400;  % Average Earth-Moon Distance (km)
n = 29.5306 * 86400; % Synodic Earth-Moon Period (s)

thrust = u * D * n^2 * m;

end