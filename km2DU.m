function DU = km2DU(km)
% I: km - Distance in kilometers (can be scalar or array)
% O: DU - Distance in DU (nondimensionalized using Earth-Moon distance)

R_EM = 384400;  % Average Earth-Moon Distance (km)

% Convert to DU
DU = km ./ R_EM;

end