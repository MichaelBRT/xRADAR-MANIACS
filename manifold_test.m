
clear, clc, close all

% Initial parameters (same as in main, this can be trimmed for efficiency)
% plot colours
tic;


blue = [0.07,0.62,1.00];
orange = [0.988, 0.38, 0];
purple = [0.659, 0, 1];

global mu
mu = 1.215058560962404E-2;


C = 3.15;

orbit1_file = fullfile('filtered PlanarOrbitData/(3,2) Cycler.csv');
%orbit2_file = "PlanarOrbitData\Lyapunov (L2).csv";

% Extracting arrays of orbit data from JPL database
i_data = readmatrix(orbit1_file);
%f_data = readmatrix(orbit2_file);

% Finding closest Jacobi Constant
ind_i = findClosestJacobi(i_data,C);
%ind_f = findClosestJacobi(f_data,C);


% Extract initials state of orbits in question
X0i = i_data(ind_i, [2,3,5,6])';
%X0f = f_data(ind_f, [2,3,5,6])';

%
Ti = i_data(ind_i, 9);
%Tf = f_data(ind_f, 9);

% ODE function ()
cr3bp = @(t,x) CR3BPMC2D(x,mu);
cr3bp_opts = odeset('RelTol',3e-13,'AbsTol',1e-13);

[~,Xi] = ode45(cr3bp,[0,Ti],X0i,cr3bp_opts);
xi = Xi(:,1); yi = Xi(:,2);
%[~,Xf] = ode45(cr3bp,[0,Tf],X0f,cr3bp_opts);
%xf = Xf(:,1); yf = Xf(:,2);

N = 60;
[Xu,Wup] = get_unst_manifold(orbit1_file,C,N,1);
%[Xs,Wsm] = get_unst_manifold(orbit2_file,C,N,-1,-1);

for i = 1:length(Xu)
    Xsi = Xu{i};
    Xsi(:,2:3) = -Xsi(:,2:3);
    Xs{i} = Xsi;
end
Wsm = Wup;
Wsm(:,2:3) = -Wup(:,2:3);


f1 = figure();
f1_ax = gca;
hold (gca, 'on')


for i = 1:N
    X = Xu{i};
    plot(f1_ax, X(:,1),X(:,2),'r')
    hold on
    X = Xs{i};
    plot(f1_ax, X(:,1),X(:,2),'g')
end

plot(xi,yi,'LineWidth',2,'Color',blue)
%plot(xf,yf,'LineWidth',2,'Color',purple);


drawEarthMoonSystem(f1_ax,1,i_data(ind_i,8));

f1.Theme = 'dark';

xlim([0.759 1.200])
ylim([-0.169 0.179])


% Search for "closest" crossing
diffs_Mat = zeros(size(Wup,1),size(Wsm,1));

for i = 1:size(Wup,1)
    for j = 1:size(Wsm,1)

        diffs_Mat(i,j) = sqrt((Wup(i,2)-Wsm(j,2))^2 +...
            0.01*(Wup(i,4)-Wsm(j,4))^2 + 0.001*(Wup(i,3)-Wsm(j,3))^2);

    end
end
%
[~, linearIdx] = min(diffs_Mat(:));
[rowIdx, colIdx] = ind2sub(size(diffs_Mat), linearIdx);

% rowIdx - index of crossing on the unstable manifold
% colIdx - index of crossing on the stable manifold

f2 = figure(2);
f2.Theme = 'dark';
hold on

scatter(Wup(:,1),Wup(:,3),7.5,'r','filled')
plot(Wup(rowIdx,1),Wup(rowIdx,4),'Color','w','Marker','x','LineWidth',1.5 ...
    ,'MarkerSize',10)
plot(Wsm(colIdx,1),Wsm(colIdx,3),'Color','w','Marker','x','LineWidth',1.5 ...
    ,'MarkerSize',10)
scatter(Wsm(:,1),Wsm(:,3),7.5,'g','filled')

grid on


toc; % first recording - how long does manifold generation take?

function [X,PHI] = stm(tau,tspan,Yspan)
% Interpolates state transition matrix from propagation along single period
% of a PO.

T = tspan(end);
t = tau*T;
Y = interp1(tspan,Yspan(:,1:16),t);
X = interp1(tspan,Yspan(:,17:20),t)';
PHI = reshape(Y,[4,4]);

end

function [all_Wu, Wu] = get_unst_manifold(orbit_file,C,N,sgn2)
global mu
if isempty(gcp("nocreate"))
    parpool(14);
end
vareqn = @(t,x) var2D(t,x,mu);
varopt = odeset('RelTol',3e-10,'AbsTol',1e-10);

cr3bp = @(t,x) CR3BPMC2D(x,mu);
cr3bp_opts = odeset('RelTol',3e-10,'AbsTol',1e-10, ...
    'Events',@(t,x) eX(t,x,mu));

ep = 1e-5;
T_int = 25;

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
parfor i = 1:N
    [~,all_Wu{i}, te, We,ie] = ode45(cr3bp ...
        , [0,T_int],Wu_0(:,i),cr3bp_opts);

    if ie == 1 
        Wu = [Wu; We];
    end

end


end

function [all_Ws, Ws] = get_stab_manifold(orbit_file,C,N, sgn)
global mu
if isempty(gcp("nocreate"))
    parpool(14);
end
vareqn = @(t,x) var2D(t,x,mu);
varopt = odeset('RelTol',3e-10,'AbsTol',1e-10);


cr3bp = @(t,x) CR3BPMC2D(x,mu);
cr3bp_opts = odeset('RelTol',3e-13,'AbsTol',1e-13, ...
    'Events',@(t,x) Poincare1(t,x,mu,sgn));

ep = 1e-5;


tau = linspace(0,1,N);


data = readmatrix(orbit_file);
idx = findClosestJacobi(data,C);
X0 = data(idx, [2,3,5,6])';
T = data(idx,9);
phi_0 = reshape(eye(4),[16,1]);

T_int = 20;

Y0 = [phi_0;X0];

[~,Yspan] = ode113(vareqn,T*tau,Y0,varopt);
M = reshape(Yspan(end,1:16)',[4,4]);

[V,D] = eigs(M);


all_eig_shifted = diag(D);

[~,important_idx] = maxk(all_eig_shifted,2,'ComparisonMethod','abs');

ind = important_idx(2);
u = real(V(:,ind));

Ws_0 = zeros(4,N);

for i = 1 :length(tau)
    X = Yspan(i,17:20)';
    PHI = reshape(Yspan(i,1:16)',[4,4]);


    v = PHI*u;
    v = u/norm(u);

    %Wsp_0 = X + ep*u;
    Ws_0(:,i) = X - ep*v;
end

Ws = [];
parfor i = 1:N
    [~,all_Ws{i}, te, We,ie] = ode45(cr3bp ...
        , [0,-T_int],Ws_0(:,i),cr3bp_opts);
    if abs(te) < 10
        if ie == 1
            Ws = [Ws; We];
        end
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


if  abs(x(4)) < 5
    value = x(1)-(1-mu);
    isterminal = 0;

else
    value = [];
    isterminal=0;
end


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
x = linspace(-1.5,1.5,1000);
y = linspace(-1.5,1.5,1000);
[X,Y] = meshgrid(x,y);

Omega = @(x,y) 1/2*(x.^2 + y.^2) + (1-mu)./sqrt((x+mu).^2 + y.^2)+...
    mu./sqrt((x-1+mu).^2 + y.^2);
hold(ax,"on")
contour(ax, X,Y,Omega(X,Y),[C/2,C/2],'w','LineWidth',2)
hold(ax, "off")
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