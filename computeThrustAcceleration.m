% To Compute thrust acceleration vector using CR3BP dynamics
% inputs 
% All states are N*1 vectors
% outputs
% u will be N*3 matrix, the columns are ux, uy and uz, rows are different
% time steps
% thrust is n*1 vector (stores thrust mag at each time step)
function [u, thrust] = ThrustAcceleration(x, y, z, dx, dy, dz, ddx, ddy, ddz, mu)
    % 
    % Number of points
    N = length(x);
    % Preallocate
    u = zeros(N, 3);
    thrust = zeros(N, 1);
    %
    for i = 1:N
        % Current position
        xi = x(i); yi = y(i); zi = z(i);
        r1 = sqrt((xi + mu)^2 + yi^2 + zi^2);
        r2 = sqrt((xi - 1 + mu)^2 + yi^2 + zi^2);
        dUdx = xi - (1 - mu)*(xi + mu)/r1^3 - mu*(xi - 1 + mu)/r2^3;
        dUdy = yi - (1 - mu)*yi/r1^3 - mu*yi/r2^3;
        dUdz =     - (1 - mu)*zi/r1^3 - mu*zi/r2^3;
        ddxi = ddx(i); ddyi = ddy(i); ddzi = ddz(i);
        dxi = dx(i); dyi = dy(i);
        
        % Compute thrust components
        ux = ddxi - 2*dyi - dUdx;
        uy = ddyi + 2*dxi - dUdy;
        uz = ddzi - dUdz;
        u(i,:) = [ux, uy, uz];
        thrust(i) = norm([ux, uy, uz]);
    end
end
