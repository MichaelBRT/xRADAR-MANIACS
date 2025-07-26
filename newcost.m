function J = newcost(X, params)
global n scalefactor
s=X(end);
Z=X(1:end-1);
check = params.is_T_design;
M = params.M;
if check
    T = exp(Z(end))/scalefactor;
else
    T = params.T;
end


tau = linspace(0,1,M);

td = tau*T/n;
a = getAccel2D(Z,params);

rho = 100;
J = trapz(td,a) + rho*s;
end