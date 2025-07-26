clear, clc, close all


rho = 0.01215;


%% Input data
% halo orbit initial condition
x0_nominal = [7.1635002736262687e-1; 0; 0; 0; 6.6125823335399858e-1; 0];  

% halo orbit period
Tperiod    = 4.6664370637490915;


% The time span of interest is
T = 1.4; % this seems arbitrary.
tspan = [0,T];


% Kolmanovsky defines initial condition
x0 = x0_nominal + [-0.05;-0.05;0.05;0;0;0];


% He defines an ode function for easy use
u_interp = @(t) [0; 0; 0];                                          % no control input
odefun = @(t, x) cr3bp_controlled_dynamics(t, x, rho, u_interp );


% He then propagates the function forward to get the target position:
[t_sim, x_sim] = ode45(odefun, tspan, x0_nominal);
xf = x_sim(end,:)';


odefun_adj = @(t,x) optimalTransferV1(t,x,rho);
bcfun = @(Xa, Xb) bcoptimalTransferV1(Xa,Xb, xf,x0);

% To use the bvp4c solver, we need to initialize and define the proper
% functions.






%% Functions


function dxdt = cr3bp_controlled_dynamics(t, x, mu, u_func)
    r1 = norm([x(1)+mu; x(2); x(3)]);
    r2 = norm([x(1)-1+mu; x(2); x(3)]);

    Ux = x(1) - (1-mu)*(x(1)+mu)/r1^3 - mu*(x(1)-1+mu)/r2^3;
    Uy = x(2) - (1-mu)*x(2)/r1^3 - mu*x(2)/r2^3;
    Uz = -(1-mu)*x(3)/r1^3 - mu*x(3)/r2^3;

    u = u_func(t);
    dxdt = zeros(6,1);
    dxdt(1:3) = x(4:6);
    dxdt(4) = 2*x(5) + Ux + u(1);
    dxdt(5) = -2*x(4) + Uy + u(2);
    dxdt(6) = Uz + u(3);
end


function [Xdot]=optimalTransferV1(t,X,rho)


% States
x = X(1:6);  % [x, y, z, vx, vy, vz]

% Adjoint variables
p  = X(7:12); % adjoint variables

% Dynamics
r1 = norm([x(1) + rho, x(2), x(3)]);
r2 = norm([x(1) - 1 + rho, x(2), x(3)]);

Ux = x(1) - (1-rho)*(x(1)+rho)/r1^3 - rho*(x(1)-1+rho)/r2^3;
Uy = x(2) - (1-rho)*x(2)/r1^3 - rho*x(2)/r2^3;
Uz = -(1-rho)*x(3)/r1^3 - rho*x(3)/r2^3;

dxdt = [x(4);
    x(5);
    x(6);
    2*x(5) + Ux;
    -2*x(4) + Uy;
    Uz];

% Optimal control: u = -p(4:6)
u         = -p(4:6);
dxdt(4:6) = dxdt(4:6) + u;

dHdx = compute_dHdx(x,p, rho);
Xdot = [dxdt; -dHdx];
Xdot      = Xdot(:);
end
     
function bc = bcoptimalTransferV1(XA,XB,xT,x0)
% used with bvp4c
% initial conditions
bc(1:6) = x0(:) - XA(1:6);

% final conditions
bc(7:12) = XB(1:6)-xT(:);

end

function dHdx = compute_dHdx(x,p, rho)
x1=x(1);x2=x(2);x3=x(3);x4=x(4);x5=x(5);x6=x(6);
p1=p(1);p2=p(2);p3=p(3);p4=p(4);p5=p(5);p6=p(6);

dHdx =[(2*p4*rho - p4 - p4*x1 - 3*p5*x2 - 3*p6*x3 - p4*rho^2 + 2*p4*x1^2 - p4*x2^2 - p4*x3^2 + p4*((rho + x1 - 1)^2 + x2^2 + x3^2)^(5/2) + p4*rho*x1 + 3*p5*rho*x2 + 3*p6*rho*x3 + 3*p5*x1*x2 + 3*p6*x1*x3)/(((rho + x1 - 1)^2 + x2^2 + x3^2)^(1/2)*(rho^2 + 2*rho*x1 - 2*rho + x1^2 - 2*x1 + x2^2 + x3^2 + 1)^2);
       (2*p5*rho - p5 + 2*p5*x1 - p5*rho^2 - p5*x1^2 + 2*p5*x2^2 - p5*x3^2 + p5*((rho + x1 - 1)^2 + x2^2 + x3^2)^(5/2) - 2*p5*rho*x1 + 3*p4*x1*x2 + 3*p6*x2*x3)/(((rho + x1 - 1)^2 + x2^2 + x3^2)^(1/2)*(rho^2 + 2*rho*x1 - 2*rho + x1^2 - 2*x1 + x2^2 + x3^2 + 1)^2);
       -(p6*rho^2 + 2*p6*rho*x1 - 2*p6*rho + p6*x1^2 - 3*p4*x1*x3 - 2*p6*x1 + p6*x2^2 - 3*p5*x2*x3 - 2*p6*x3^2 + p6)/(((rho + x1 - 1)^2 + x2^2 + x3^2)^(1/2)*(rho^2 + 2*rho*x1 - 2*rho + x1^2 - 2*x1 + x2^2 + x3^2 + 1)^2);
        p1 - 2*p5;
        p2 + 2*p4;
        p3];

end
