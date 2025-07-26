clear, clc, close all

% parameters definition

mu = 0.0122;
tol = 1e-10;

odeopts = odeset('RelTol',3*tol,'AbsTol',tol);

% States from Chase and Saevar

X3 = [1.03023361221087; 
    0.0687372443157935; 
    0.147641949569388; 
    0.157623720622668; 
    -0.181091773654881; 
    0.389575654309805];

X4 = [-0.884299373220872; 
    -0.434399036504677; 
    0.337924353724463; 
    -0.382566367054743; 
    0.0896125353559663; 
    0.302512913924565];

T = 3.540885;


dynamics = @(t,x) cr3bp(t,x,mu);

[~,x1] = ode45(dynamics,[0,2*T],X3,odeopts);
[~,x2] = ode45(dynamics,[0,-T],X4,odeopts);


f = figure();
ax = gca;
f.Theme = 'dark';

hold(ax,'on')
plot3(ax,x1(:,1),x1(:,2),x1(:,3),'LineWidth',2)
plot3(ax,x2(:,1),x2(:,2),x2(:,3),'LineWidth',2)
drawEarthMoonSystem3D(mu,ax,1);
hold(ax,'off')


f2 = figure();
ax2 = gca;
f2.Theme = 'dark';

hold(ax,'on')
plot(ax2,x1(:,1),x1(:,2),'LineWidth',2)
plot(ax2,x2(:,1),x2(:,2),'LineWidth',2)
drawEarthMoonSystem(ax,1,mu,C)
hold(ax,'off')







%% 3-Dimensional CR3BP Dynamics
function C = Jconst(X,mu)
    x = X(1);
    y = X(2);
    z = X(3);

    vx = X(4);
    vy = X(5);
    vz = X(6);

    r1 = sqrt((x+mu)^2 + y^2 + z^2);
    r2 = sqrt((x-1+mu)^2 + y^2 + z^2);

    Ubar = -1/2*(x^2+y^2+z^2) - (1-mu)/r1 -mu/r2;

    C = -(vx^2 + vy^2 + vz^2) -2*Ubar;

end


function xd = cr3bp(~,X,mu)
    xd = zeros(6,1);
    
    x = X(1);
    y = X(2);
    z = X(3);

    vx = X(4);
    vy = X(5);
    vz = X(6);

    r1 = sqrt((x+mu)^2 + y^2 + z^2);
    r2 = sqrt((x-1+mu)^2 + y^2 + z^2);

    xd(1) = vx;
    xd(2) = vy;
    xd(3) = vz;

    xd(4) = 2*vy + x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3;
    xd(5) = -2*vx + y - (1-mu)*y/r1^3 - mu*y/r2^3;
    xd(6) = -(1-mu)*z/r1^3 - mu*z/r2^3;

end


function drawEarthMoonSystem(ax,mode,mu,C)
% mode: 0 - light mode
%       1 - dark mode


blue = [0.07,0.62,1.00];

if nargin < 2
    mode = 0;
end

if mode == 1
    moon_colour = [0.8,0.8,0.8];
else
    moon_colour = [0.2,0.2,0.2];
end

% Radii to scale
Re = 6378/389703;
Rm = 1737/389703;

% Lagrange point positions
xL1 = Lptpos(mu,1);
xL2 = Lptpos(mu,2);
xL3 = Lptpos(mu,3);

%Jconst([xL3;0;0;0])

hold(ax,"on")
utils.draw_shaded_circle(ax,[mu,0],Re, blue,1)
utils.draw_shaded_circle(ax,[1-mu,0],Rm, moon_colour,1)
plot(ax,xL1, 0, 'rx')
plot(ax,xL2, 0, 'rx')
plot(ax,xL3, 0, 'rx')
plot(ax,1/2,sqrt(3)/2,'rx','LineWidth',1)
plot(ax,1/2,-sqrt(3)/2,'rx','LineWidth',1)
if nargin > 2
    ZVC_no_shade(C,mu,ax,mode)
end

hold (ax,'off')
xlabel(ax,'$x$ [DU]','Interpreter','latex')
ylabel(ax,'$y$ [DU]','Interpreter','latex')
set(ax,'fontsize',16)

grid on
axis equal
end


function ZVC_no_shade(C,mu,ax,mode)
x = linspace(-1.5,1.5,1000);
y = linspace(-1.5,1.5,1000);
[X,Y] = meshgrid(x,y);

Omega = @(x,y) 1/2*(x.^2 + y.^2) + (1-mu)./sqrt((x+mu).^2 + y.^2)+...
    mu./sqrt((x-1+mu).^2 + y.^2);
hold(ax,"on")
if mode == 1
contour(ax, X,Y,Omega(X,Y),[C/2,C/2],'w','LineWidth',2)
else
contour(ax, X,Y,Omega(X,Y),[C/2,C/2],'k','LineWidth',2)
end
hold(ax, "off")
end


function drawEarthMoonSystem3D(mu,ax,mode)
if mode == 1
    moon_color = [0.9,0.9,0.9];
else
    moon_color = [0.1,0.1,0.1];
end


%% plotting Earth, Moon and Lagrange Points
Re = 6378/389703;
Rm = 1737/389703;


[Xe,Ye,Ze] = ellipsoid(ax,-mu,0,0,Re,Re,Re,400);
[Xm,Ym,Zm] = ellipsoid(ax,1-mu,0,0,Rm,Rm,Rm,400);

xL1 = Lptpos(mu,1);
xL2 = Lptpos(mu,2);
xL3 = Lptpos(mu,3);

hold (ax, 'on')
surf(ax,Xe,Ye,Ze,'EdgeColor','b','FaceColor','flat')
surf(ax,Xm,Ym,Zm,'EdgeColor',moon_color,'FaceColor','flat');
plot3(ax, xL1, 0, 0, 'rx')
plot3(ax, xL2, 0, 0, 'rx')
plot3(ax, xL3, 0, 0, 'rx')
plot3(ax,1/2,sqrt(3)/2,0,'rx','LineWidth',1)
plot3(ax,1/2,-sqrt(3)/2,0,'rx','LineWidth',1)
hold (ax,'off')


xlabel(ax,'$x$ [DU]','FontSize',14,'Interpreter','latex')
ylabel(ax,'$y$ [DU]','FontSize',14,'Interpreter','latex')
zlabel(ax,'$z$ [DU]','FontSize',14,'Interpreter','latex')
title(ax,'')
axis(ax,'equal')
view(ax,[0,0,1])
grid(ax,"on")
box(ax,"on")
pbaspect(ax, [1,1,1]);

end
