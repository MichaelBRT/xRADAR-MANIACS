clear, clc, close all
blue     = [0.07, 0.62, 1.00];
orange   = [0.988, 0.38, 0];
purple   = [0.659, 0, 1];
gray     = [0.1, 0.1, 0.1];





load('scenarios\32Res-L1Lya_eX.mat')

X1 = Xi;

f = figure();
hold on

plot(init_arc(:,1),init_arc(:,2),'-','LineWidth',2,'Color',blue)
dVt = 1000*deltaV_req;

load('scenarios\L1Lya_to_DPO.mat')
plot(init_arc(:,1),init_arc(:,2),'-','LineWidth',2,'Color',[0,0.5,0])
dVt = dVt + 1000*deltaV_req;


load('scenarios\DPO_to_31res.mat')
plot(init_arc(:,1),init_arc(:,2),'-','LineWidth',2,'Color',orange)
dVt = dVt + 1000*deltaV_req;
X2 = Xf;

plot(X1(:,1),X1(:,2),'LineWidth',2,'Color','k')
plot(X2(:,1),X2(:,2),'LineWidth',2,'Color','k')
drawEarthMoonSystem(gca(),0);
title('Earth-Moon Synodic Frame','FontSize',16,'Interpreter','latex')

% legend(initOrbitName,targetOrbitName,'Transfer Trajectory', ...
%     'Thrusting Phase','Thrust vectors')
% 




%% Functions

function [all_Wu, Wu,tu] = get_unst_manifold(orbit_file,C,N,sgn2)
global mu
if isempty(gcp("nocreate"))
    parpool(14);
end
vareqn = @(t,x) var2D(t,x,mu);
varopt = odeset('RelTol',3e-10,'AbsTol',1e-10);

cr3bp = @(t,x) CR3BPMC2D(x,mu);
cr3bp_opts = odeset('RelTol',3e-10,'AbsTol',1e-10, ...
    'Events',@(t,x) mY(t,x,mu));

ep = 1e-5;
T_int = 10;

tau = linspace(0,1,N);


data = readmatrix(orbit_file);
idx = findClosestJacobi(data,C);
X0 = data(idx, [2,3,5,6])';
T = data(idx,9);

phi_0 = reshape(eye(4),[16,1]);

Y0 = [phi_0;X0];

[~,Yspan] = ode113(vareqn,T*tau,Y0,varopt);


M = reshape(Yspan(end,1:16)',[4,4]);
[V,D] = eigs(M);

all_eig_shifted = diag(D);

[~,important_idx] = maxk(all_eig_shifted,2,'ComparisonMethod','abs');

ind = important_idx(1);
u = real(V(:,ind));


%Wsp_0 = zeros(4,N);
Wu_0 = zeros(4,N);



for i = 1 :length(tau)
    X = Yspan(i,17:20)';
    PHI = reshape(Yspan(i,1:16)',[4,4]);

    v = PHI*u;
    v = sgn2*v/norm(v);

    Wu_0(:,i) = X + ep*v;
end

Wu = [];
tu = [];
parfor i = 1:N
    [~,all_Wu{i}, te, We,ie] = ode45(cr3bp ...
        , [0,T_int],Wu_0(:,i),cr3bp_opts);

    if ie == 1
        Wu = [Wu; We];
        tu = [tu;te];
    end

end


end



function [value, isterminal, direction] = eY(t,x,mu)
% v = [x(3);x(4)];
% phi = acosd(v'*[1;0]);

% r2 = sqrt((x(1)-1+mu)^2 + x(2)^2);
% Rm = 1737/389703;


if  abs(x(4)) < 5
    value = x(1)-(-mu);
    isterminal = 0;

else
    value = [];
    isterminal=0;
end


end


function [value, isterminal, direction] = mY(t,x,mu)
% v = [x(3);x(4)];
% phi = acosd(v'*[1;0]);

% r2 = sqrt((x(1)-1+mu)^2 + x(2)^2);
% Rm = 1737/389703;


%if  abs(x(4)) < 5
value = x(1)-(1-mu);
isterminal = 1;

% else
%     value = [];
%     isterminal=0;
% end
direction = 0;


end



function [value, isterminal, direction] = eX(t,x,mu)
% v = [x(3);x(4)];
% phi = acosd(v'*[1;0]);

% r2 = sqrt((x(1)-1+mu)^2 + x(2)^2);
% Rm = 1737/389703;

xL1 = Lptpos(mu,1);

if  x(1) < xL1 && x(1) > -1.5 && abs(x(3)) < 5
    value = x(2);
    isterminal = 0;

else
    value = [];
    isterminal=0;
end

direction = 0;

end

function [value, isterminal, direction] = mX(t,x,mu)
% v = [x(3);x(4)];
% phi = acosd(v'*[1;0]);

% r2 = sqrt((x(1)-1+mu)^2 + x(2)^2);
% Rm = 1737/389703;

xL1 = Lptpos(mu,1);

if  x(1) > xL1 && abs(x(1)) < 3 && abs(x(3)) < 5
    value = x(2);
    isterminal = 0;

else
    value = [];
    isterminal=0;
end




direction(1) = 0;

end






%% Plot Zero Velocity Curves
function ZVC(C,mu,ax)
%ax = app.ax;
x = linspace(-1.5,1.5,1000);
y = linspace(-1.5,1.5,1000);
[X,Y] = meshgrid(x,y);

Omega = @(x,y) 1/2*(x.^2 + y.^2) + (1-mu)./sqrt((x+mu).^2 + y.^2)+...
    mu./sqrt((x-1+mu).^2 + y.^2);

Z = Omega(X,Y);

hold(ax, 'on');

% ZVC contour boundary at C/2
Clevel = C / 2;

contourData = contourc(x, y, Z, [Clevel Clevel]);

% Extract contour polygon(s)
idx = 1;
while idx < size(contourData, 2)
    numPoints = contourData(2, idx);
    xPoly = contourData(1, idx+1:idx+numPoints);
    yPoly = contourData(2, idx+1:idx+numPoints);

    % Fill each forbidden region
    patch(ax, xPoly, yPoly, [0.5,0.5,0.5], ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
        'DisplayName', 'ZVC Forbidden', ...
        'Tag', 'ZVC');
    idx = idx + numPoints + 1;
end


% Draw the ZVC line on top
contour(ax, X, Y, Z, [Clevel Clevel], ...
    'LineColor', 'k', 'LineWidth', 0.35, ...
    'Tag', 'ZVC');
end


%% Find Nearest Jacobi Constant
function idx = findClosestJacobi(data, targetC)
[~, idx] = min(abs(data(:,8) - targetC));
end

%% Draw shaded circle of specified size and color
function draw_shaded_circle(ax, center, radius, color, fillCircle)
if nargin < 4, color = [0.3, 0.6, 1]; fillCircle = true; end
theta = linspace(0, 2*pi, 50);
xCirc = center(1) + radius * cos(theta);
yCirc = center(2) + radius * sin(theta);
if fillCircle
    fill(ax, xCirc, yCirc, color, 'EdgeColor', 'none');
else
    plot(ax, xCirc, yCirc, 'Color', color, 'LineStyle', '-.', 'LineWidth', 0.5);
end
end

%% Variational E
function PHIdot = var2D(~,PHI, mu)
% PHIdot=var2D(t,PHI)
%
% This here is a preliminary state transition, PHI(t,t0),
% matrix equation attempt for the planar CR3BP, based on...
%
%        d PHI(t, t0)
%        ------------ =  F(t) * PHI(t, t0)
%             dt
%
%-----------------------------------------------------------
% CR3BP CONVENTION:
%                 L4
%
%
%    L3-----M1-------L1---M2---L2         M1=1-mu, M2=mu
%
%
%                 L5
%
% Shane Ross (revised 9.23.97)
%
% global FORWARD mu


mu1=1-mu;
mu2=  mu;

x(1:4) = PHI(17:20);
phi  = reshape(PHI(1:16),4,4);

r2= (x(1)+mu )^2 + x(2)^2;	% r: distance to m1, LARGER MASS
R2= (x(1)-mu1)^2 + x(2)^2;	% R: distance to m2, smaller mass
r3= r2^1.5; r5= r2^2.5;
R3= R2^1.5; R5= R2^2.5;

omgxx= 1+(mu1/r5)*(3*(x(1)+mu2)^2)+(mu2/R5)*(3*(x(1)-mu1)^2)-(mu1/r3+mu2/R3);
omgyy= 1+(mu1/r5)*(3* x(2)^2     )+(mu2/R5)*(3* x(2)^2     )-(mu1/r3+mu2/R3);
omgxy= 3*x(2)*     (mu1*(x(1)+mu2)/r5+mu2*(x(1)-mu1)/R5);

F     =[   0     0     1     0  ;
    0     0     0     1  ;
    omgxx omgxy   0     2  ;
    omgxy omgyy  -2     0 ];

phidot = F * phi; % variational equation

PHIdot        = zeros(20,1);
PHIdot(1:16)  = reshape(phidot,16,1);
PHIdot(17)    = x(3);
PHIdot(18)    = x(4);
PHIdot(19)    = x(1)-(mu1*(x(1)+mu2)/r3) -(mu2*(x(1)-mu1)/R3) + 2*x(4);
PHIdot(20)    = x(2)-(mu1* x(2)     /r3) -(mu2* x(2)     /R3) - 2*x(3);
PHIdot(17:20) = PHIdot(17:20);

end


function drawEarthMoonSystem(ax,mode,C)
% mode: 0 - light mode
%       1 - dark mode
mu = 1.215058560962404E-2;

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
draw_shaded_circle(ax,[mu,0],Re, blue,1)
draw_shaded_circle(ax,[1-mu,0],Rm, moon_colour,1)
plot(ax,xL1, 0, 'rx')
plot(ax,xL2, 0, 'rx')
plot(ax,xL3, 0, 'rx')
plot(ax,1/2,sqrt(3)/2,'rx','LineWidth',1)
plot(ax,1/2,-sqrt(3)/2,'rx','LineWidth',1)
if nargin > 2
    ZVC(C,mu,ax)
end

hold (ax,'off')
xlabel(ax,'$x$ [DU]','Interpreter','latex')
ylabel(ax,'$y$ [DU]','Interpreter','latex')
set(ax,'fontsize',16)

grid on
axis equal
end


function C = Jconst(X)
% Function calculates Jacobi constant given a state
x = X(1);
y = X(2);
xdot = X(3);
ydot = X(4);

mu = 0.012150584270572;
mu1=1-mu;
mu2=  mu;

r1 = ((x+mu2)^2 + y^2)^(1/2);
r2 = ((x-mu1)^2 + y^2)^(1/2);

U = -1/2*(x^2+y^2)-mu1/r1-mu2/r2; %-1/2*mu1*mu2;

C = -(xdot^2 + ydot^2)-2*U;
end

function [idx_u, idx_s] = get_manifold_intersection(Wu,Ws)
diffs_Mat = zeros(size(Wu,1),size(Ws,1));

a = 0.1;
b = 1;
c = 1;
for i = 1:size(Wu,1)
    for j = 1:size(Ws,1)
        diffs_Mat(i,j) = sqrt(a*(Wu(i,2)-Ws(j,2))^2 +...
            b*(Wu(i,4)-Ws(j,4))^2 + c*(Wu(i,3)-Ws(j,3))^2);


    end
end

[~, linearIdx] = min(diffs_Mat(:));
[idx_u, idx_s] = ind2sub(size(diffs_Mat), linearIdx);
end