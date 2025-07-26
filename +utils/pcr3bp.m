function xd = pcr3bp(~,X,mu)
    xd = zeros(4,1);
    
    x = X(1);
    y = X(2);
    vx = X(3);
    vy = X(4);

    r1 = sqrt((x+mu)^2 + y^2);
    r2 = sqrt((x-1+mu)^2 + y^2);

    xd(1) = vx;
    xd(2) = vy;
    xd(3) = 2*vy + x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3;
    xd(4) = -2*vx + y - (1-mu)*y/r1^3 - mu*y/r2^3;

end