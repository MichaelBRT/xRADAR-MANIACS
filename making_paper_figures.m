clear, clc, close all
mu =  0.01215058560962404;


f = figure();
ax = gca;

C = 3.15;
utils.drawEarthMoonSystem(ax,0)
grid(ax,'off')

%xax = xline(ax,0,'LineWidth',1,'Color','k');
%yax = yline(ax,0,'LineWidth',1,'Color','k');
ZVC(C,mu,ax)
%ax.Children = flip(ax.Children);
drawnow;
axis(ax,'equal')
axis(ax,'off')


%%
clear, clc, close all
global mu
mu =  0.01215058560962404;
C = 3.165;

blue     = [0.07, 0.62, 1.00];
orange   = [0.988, 0.38, 0];
purple   = [0.659, 0, 1];
gray     = [0.1, 0.1, 0.1];

orbit1_file = fullfile('filtered PlanarOrbitData/Lyapunov (L1).csv');
orbit2_file = "PlanarOrbitData\Lyapunov (L2).csv";

% Extracting arrays of orbit data from JPL database
i_data = readmatrix(orbit1_file);
f_data = readmatrix(orbit2_file);

% Finding closest Jacobi Constant
ind_i = findClosestJacobi(i_data,C);
ind_f = findClosestJacobi(f_data,C);


% Extract initials state of orbits in question
X0i = i_data(ind_i, [2,3,5,6])';
X0f = f_data(ind_f, [2,3,5,6])';

%
Ti = i_data(ind_i, 9);
Tf = f_data(ind_f, 9);

% ODE function ()
cr3bp = @(t,x) CR3BPMC2D(x,mu);
cr3bp_opts = odeset('RelTol',3e-13,'AbsTol',1e-13);

[~,Xi] = ode45(cr3bp,[0,Ti],X0i,cr3bp_opts);
xi = Xi(:,1); yi = Xi(:,2);
[~,Xf] = ode45(cr3bp,[0,Tf],X0f,cr3bp_opts);
xf = Xf(:,1); yf = Xf(:,2);

N = 60;
[Xu,Wup,tu] = get_unst_manifold(orbit1_file,C,N,1);
[Xs,Wsm,ts] = get_unst_manifold(orbit2_file,C,N,-1);
Wsm(:,2:3) = -Wsm(:,2:3);
ts = -ts;



f2 = figure();
f2_ax = gca;
hold (gca, 'on')


[iU,iS] = get_manifold_intersection(Wup,Wsm);


for i = 1:N
    X = Xu{i};
    plot(f2_ax, X(:,1),X(:,2),'Color',orange,'LineWidth',1)
    X = Xs{i};
    plot(f2_ax, X(:,1),-X(:,2),'Color',purple,'LineWidth',1)
    hold on
end
X = Xu{iU};

plot(f2_ax, X(:,1),X(:,2),'Color',blue,'LineWidth',3)
X = Xs{iS};
plot(f2_ax, X(:,1),-X(:,2),'Color',[0,0.5,0],'LineWidth',3)
X = Xu{iU};


axis(f2_ax,'off')

plot(xi,yi,'LineWidth',2,'Color','k')
plot(xf,yf,'LineWidth',2,'Color','k')
plot(xi(1),yi(1),'o','LineWidth',2,'Color','k','MarkerSize',6 ...
    ,'MarkerFaceColor','w')
plot(xf(1),yf(1),'o','LineWidth',2,'Color','k','MarkerSize',6 ...
    ,'MarkerFaceColor','w')

plot(f2_ax, X(1,1),X(1,2),'o','LineWidth',3,'Color',blue,'MarkerSize',8 ...
    ,'MarkerFaceColor','w')
X = Xs{iS};

plot(f2_ax, X(1,1),-X(1,2),'o','LineWidth',3,'Color',[0,0.5,0],'MarkerSize',8 ...
    ,'MarkerFaceColor','w')

utils.draw_shaded_circle(f2_ax,[1-mu,0],1737/389703,'k',1)
axis(f2_ax,'equal')


f3 = figure();
f3_ax = gca();
hold(f3_ax,'on')
plot(f3_ax,Wup(:,2),Wup(:,4),'Color', orange,'LineWidth',2')
plot(f3_ax,Wsm(:,2),Wsm(:,4),'Color', purple,'LineWidth',2')
plot(f3_ax,Wup(iU,2),Wup(iU,4),'Color', blue,'LineWidth',2','Marker','x' ...
    ,'MarkerSize',15)
plot(f3_ax,Wsm(iS,2),Wsm(iS,4),'Color', [0,0.5,0],'LineWidth',2','Marker','x' ...
    ,'MarkerSize',15)
grid(f3_ax,"on")
xlabel('$y$ [DU]','FontSize',16,'Interpreter','latex')
ylabel('$\dot{y}$ [DU/TU]','FontSize',16,'Interpreter','latex')
f3_ax.XLim = [-0.0985   -0.0011];
f3_ax.YLim = [-1.4414    2.3866];
%f3_ax.Position = [0.1615    0.1100    0.7435    0.8150];
%axis(f2_ax,'tight')


%%
tol = 1e-10;
cr3bp = @(t,x) utils.pcr3bp(t,x,mu);

mu =  0.01215058560962404;
% Initial parameters (same as in main, this can be trimmed for efficiency)
T_em = 2.361 * 10^6;
D = 3.850*10^5;
n = 2*pi/T_em;
TU_to_Days = T_em/(2*pi*3600*24);
mass = 1000;


M = 1000;
tau = linspace(0,1,M);

rowIdx = iU;
colIdx = iS;
Wu = Wup;
Ws = Wsm;

cr3bp_opts = odeset('RelTol',3e-14,'AbsTol',1e-14);
XX0i = Wu(rowIdx,:)';
XX0f = Ws(colIdx,:)';

% timespans of ballistic arcs
TTi = tu(rowIdx);
TTf = ts(colIdx);

% span
[tti,XXi] = ode113(cr3bp,[0,-TTi],XX0i,cr3bp_opts);
[ttf,XXf] = ode113(cr3bp,[0,-TTf],XX0f,cr3bp_opts);

% Approximate time of the total transfer
TOF = abs(TTi)+abs(TTf);

% Assume thrusting arc takes up some fraction of time of the total transfer
% do a sweep of maneuver time to locate either:
% -maneuver time that yields minimum delta V
% -maneuver time that yields minimum max thrust
frac_span = linspace(0.01,0.5,50);

% Creating Maneuver Arc
for j = 1:length(frac_span)
    frac = frac_span(j);

    % Time of total transfer (ballistic + manuever) (estimate)


    % Time of just maneuver
    TOM = frac*TOF;
    TOM_span(j) = TOM;


    % Initial and final states of maneuver arc
    Si = interp1(tti,XXi,-TOM/2);
    Sf = interp1(ttf,XXf, TOM/2);


    % propagate from the crossing (forwards on unstable manifold, backwards
    % on stable manifold) the length of the maneuver
    [ti,si] = ode113(cr3bp,[0,TOM],Si,cr3bp_opts);
    [tf,sf] = ode113(cr3bp,[0,-TOM],Sf,cr3bp_opts);

    % flip the stable manifold arc
    tf = flip(tf) + TOM;
    sf = flip(sf);

    % interpolate the two arcs at M discrete points
    ssi = interp1(ti/TOM,si,tau);
    ssf = interp1(tf/TOM,sf,tau);

    % generate a third arc that smoothly moves from one arc to the other
    [S0,dS0] = utils.interp_arc(ssi,ssf,M,TOM,mu);

    % Save the intermediate arc
    all_thrust_arcs{j} = S0;

    x = S0(:,1); y = S0(:,2);
    vx = S0(:,3); vy = S0(:,4);
    ax = dS0(:,3); ay = dS0(:,4);

    r1 = sqrt((x+mu).^2 + y.^2);
    r2 = sqrt((x-1+mu).^2 + y.^2);

    ux{j} = ax - 2*vy - x + (1-mu)*(x+mu)./r1.^3 + mu*(x-1+mu)./r2.^3;
    uy{j} = ay + 2*vx - y + (1-mu)*y./r1.^3 + mu*y./r2.^3;

    u = sqrt(ux{j}.^2 + uy{j}.^2);

    % dimensionalize acceleration
    a = n^2*D*u;   % km/s^2

    % Save thrust
    all_thrust{j} = mass*a*1000; %N


    td = tau*TOM/n;

    td_cell{j} = td;

    Thr_max(j) = max(mass*a);
    dV(j) = trapz(td,a);

end

% Extract important information
% indices for either min dV or min max thrust
[~,dVmin_idx] = min(dV);
[~, Thrmin_idx] = min(Thr_max);




idx = dVmin_idx;




thrust_arc = all_thrust_arcs{idx};
thrust_time = frac_span(idx)*TOF;

unstable_idx = find(abs(tti)>=thrust_time/2);
unstable_arc = flip(XXi(unstable_idx,:));

stable_idx = find(abs(ttf)>=thrust_time/2);
stable_arc = XXf(stable_idx,:);

% entire transfer trajectory:
% unstable manifold + thrust arc + stable manifold
init_arc = [unstable_arc;thrust_arc;stable_arc];

% nondimensional normalized time corresponding to init_arc
tau = 1/TOF * linspace(0,TOF,size(init_arc, 1));

% Obtain thrust profile:
Nu = size(unstable_arc,1);
Ns = size(stable_arc,1);


ux = [zeros(Nu,1);ux{idx};zeros(Ns,1)];
uy = [zeros(Nu,1);uy{idx};zeros(Ns,1)];
u = sqrt(ux.^2 + uy.^2);
U = [ux,uy,u];

thrust = [zeros(Nu,1);all_thrust{idx};zeros(Ns,1)];

% Required delta V
deltaV_req = dV(idx);

% dimensional time corresponding to init_arc
td = (tau*TOF/n)'/(3600*24); % days
tTU = tau*TOF;



f4 = figure();
f4_ax = gca;
hold on

X = Xu{iU};

%plot(f4_ax, X(:,1),X(:,2),'Color',blue,'LineWidth',3)
X = Xs{iS};
%plot(f4_ax, X(:,1),-X(:,2),'Color',[0,0.5,0],'LineWidth',3)
X = Xu{iU};


axis(f4_ax,'off')

plot(xi,yi,'LineWidth',2,'Color','k')
plot(xf,yf,'LineWidth',2,'Color','k')
%plot(xi(1),yi(1),'o','LineWidth',2,'Color','k','MarkerSize',6 ...
%    ,'MarkerFaceColor','w')
%plot(xf(1),yf(1),'o','LineWidth',2,'Color','k','MarkerSize',6 ...
%    ,'MarkerFaceColor','w')

%plot(f4_ax, X(1,1),X(1,2),'o','LineWidth',3,'Color',blue,'MarkerSize',8 ...
%    ,'MarkerFaceColor','w')
X = Xs{iS};

%plot(f4_ax, X(1,1),-X(1,2),'o','LineWidth',3,'Color',[0,0.5,0],'MarkerSize',8 ...
%    ,'MarkerFaceColor','w')
plot([unstable_arc(:,1);thrust_arc(1,1)],[unstable_arc(:,2);thrust_arc(1,2)] ...
    ,'LineWidth',2,'Color',orange)
plot([thrust_arc(:,1);stable_arc(1,1)],[thrust_arc(:,2);stable_arc(1,2)], ...
    'LineWidth',2,'Color',blue)
plot(stable_arc(:,1),stable_arc(:,2),'LineWidth',2,'Color',purple)
utils.draw_shaded_circle(f4_ax,[1-mu,0],1737/389703,'k',1)
axis(f4_ax,'equal')


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