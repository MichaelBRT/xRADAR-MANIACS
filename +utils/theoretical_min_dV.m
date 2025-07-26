function dVmin = theoretical_min_dV(C1,C2,mu)
    DU = 3.850*10^5;
    H_LLO = 100 / DU; % low lunar orbit altitude

    Rm = 1737 / DU; % Lunar Radius
    
    x = (1-mu) + Rm + H_LLO;

    Ubar = -1/2 * x^2  - (1-mu)/abs(x+mu) - mu/abs(x-1+mu);

    V = @(C) sqrt(-2*Ubar - C);

    dVmin = 1000* abs(V(C2)- V(C1)); % m/s

end