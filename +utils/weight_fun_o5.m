function [w,dw] = weight_fun_o5(m,a,b)
    tau = linspace(0,1,m)';
    tauS = 2*tau - 1;

    phi_ab = [tauS.^7 - 3 * tauS.^3 + 2*tauS, tauS.^5 - 2*tauS.^3 + tauS];

    Dphi_ab = 2*[7*tauS.^6 - 9*tauS.^2 + 2, 5*tauS.^4-6*tauS.^2 + 1];

    phi_0 = -0.25*tauS.^3 + 0.75*tauS + 0.5;
    Dphi_0 = -0.75*tauS.^2 - + 0.75;


    w = phi_ab*[a;b] + phi_0;
    dw = Dphi_ab*[a;b] + Dphi_0;
end