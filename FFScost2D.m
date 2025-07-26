function J = FFScost2D(Z, params)
global n scalefactor scalefactor2

check = params.is_T_design;
M = params.M;
if check
    T = exp(Z(end))/scalefactor;
else
    T = params.T;
end


tau = linspace(0,1,M);

rho = params.rho;
td = tau*T/n;
a = getAccel2D(Z,params);
% s = getsval(Z,params);
s = Z(end)/scalefactor2;
J = trapz(td,a) + rho*s;
end