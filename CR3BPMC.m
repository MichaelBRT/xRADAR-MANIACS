%% Equations of Motion for CR3BP in dimensionless form
%based on slide #20 from Aerosp 548 winter 2022 Prof. Kolmanovsky
%prepared by Maxwell Chappell and Moise Mouyebe
%no external forces are currently considered

function dx = CR3BPMC(x, mu)

    dx(1:3, 1) = x(4:6); %position derivatives = velocities
    
    % rho = 1/(81.3045 + 1); %M1 = 81.3045 M2
    
    r1 = sqrt(((x(1) + mu)^2) + (x(2)^2) + (x(3)^2));
    r2 = sqrt(((x(1) - 1 + mu)^2) + (x(2)^2) + (x(3)^2));
    
    dx(4, 1) = 2*x(5) + x(1) - (((1 - mu)*(x(1) + mu))/(r1^3)) - ((mu*(x(1) - 1 + mu))/(r2^3)); %xdtdt
    
    dx(5, 1) = -2*x(4) + x(2) - (((1 - mu)*x(2))/(r1^3)) - ((mu*x(2))/(r2^3)); %ydtdt
    
    dx(6, 1) = -(((1 - mu)*x(3))/(r1^3)) - ((mu*x(3))/(r2^3)); %zdtdt

end