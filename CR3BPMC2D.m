%% Equations of Motion for CR3BP in dimensionless form
%based on slide #20 from Aerosp 548 winter 2022 Prof. Kolmanovsky
%prepared by Maxwell Chappell and Moise Mouyebe
%no external forces are currently considered

function dx = CR3BPMC2D(x, mu)

    dx(1:2, 1) = x(3:4); %position derivatives = velocities
    
    
    r1 = sqrt(((x(1) + mu)^2) + (x(2)^2));
    r2 = sqrt(((x(1) - 1 + mu)^2) + (x(2)^2));
    
    dx(3, 1) = 2*x(4) + x(1) - (((1 - mu)*(x(1) + mu))/(r1^3)) - ((mu*(x(1) - 1 + mu))/(r2^3)); %xdtdt
    
    dx(4, 1) = -2*x(3) + x(2) - (((1 - mu)*x(2))/(r1^3)) - ((mu*x(2))/(r2^3)); %ydtdt

end