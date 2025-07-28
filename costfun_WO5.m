function [J,dV,ux,uy,c] = costfun_WO5(a,b,T,Xi,Xf,mu,m)
    cr3bp = @(X) utils.pcr3bp(0,X,mu);

    tau = linspace(0,1,m);
    for i = 1:length(Xi)
        dXi(i,:) = cr3bp(Xi(i,:)')';
        dXf(i,:) = cr3bp(Xf(i,:)')';
    end

    [w,dw] = utils.weight_fun_o5(m,a,b);

    c = (1-w).*Xi + w.*Xf;
    dc = 1/T*dw.*(Xf - Xi) + (1-w).*dXi + w.*dXf;

    x = c(:,1); y = c(:,2);
    vx = c(:,3); vy = c(:,4);
    ax = dc(:,3); ay = dc(:,4);

    r1 = sqrt((x+mu).^2 + y.^2);
    r2 = sqrt((x-1+mu).^2 + y.^2);

    ux = ax - 2*vy - x + (1-mu)*(x+mu)./r1.^3 + mu*(x-1+mu)./r2.^3;
    uy = ay + 2*vx - y + (1-mu)*y./r1.^3 + mu*y./r2.^3;

    u = sqrt(ux.^2 + uy.^2);

    J = trapz(T*tau, u.^2);
    dV = trapz(T*tau, u);
      
end