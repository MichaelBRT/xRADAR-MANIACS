% Using getAmatrix and findCF to get fourier and derivatives
function [x, xd, xdd] = findFFS(tou, X, Xi, Xf, nffs)
    [A, Ad, Add] = getAmatrix(tou,nffs);
    [Cf,Cfd,Cfdd] = findCF(tou,Xi,Xf);
     x     = A * Xx + Cf;
    xd  = Ad * Xx + Cfd;
    xdd = Add * Xx + Cfdd;

    

end




