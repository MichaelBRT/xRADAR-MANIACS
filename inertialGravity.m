%% Function to get the inertial gravity acceleration terms (no vel.)
function a_g = inertialGravity(trajectory, mu)

    x = trajectory(:,1);
    y = trajectory(:,2);
    z = trajectory(:,3);
    
    r1 = sqrt(((x + mu).^2) + (y.^2) + (z.^2));
    r2 = sqrt(((x - 1 + mu).^2) + (y.^2) + (z.^2));

    a_g = zeros(length(x), 3);
    
    a_g(:, 1) = -(((1 - mu)*(x + mu))./(r1.^3)) - ((mu*(x - 1 + mu))./(r2.^3)); % a_x
    
    a_g(:, 2) = -(((1 - mu)*y)./(r1.^3)) - ((mu*y)./(r2.^3)); % a_y
    
    a_g(:, 3) = -(((1 - mu)*z)./(r1.^3)) - ((mu*z)./(r2.^3)); % a_z

end