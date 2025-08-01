global mu

mu = 1.215058560962404E-2;
% orbit_file = fullfile('filtered PlanarOrbitData/Lyapunov (L1).csv');
% data = readmatrix(orbit_file);


folder_name = 'filtered PlanarOrbitData';
folder_path = fullfile(pwd,folder_name);


files = dir(folder_path);
files = files(~[files.isdir]);  % Remove '.' and '..'

num_files = numel(files);
for k = 1:4
    switch k
        case 1
            outputfolder = fullfile('Poincaré Section Data/eX/');
            fprintf(['Current Section: Earth X \n'])
        case 2
            outputfolder = fullfile('Poincaré Section Data/mX/');
            fprintf(['Current Section: Moon X \n'])
        case 3
            outputfolder = fullfile('Poincaré Section Data/eY/');
            fprintf(['Current Section: Earth Y \n'])
        case 4
            outputfolder = fullfile('Poincaré Section Data/mY/');
            fprintf(['Current Section: Moon Y \n'])
    end
    for i = 1:num_files

        orbit_csv = files(i).name;
        orbit_file = fullfile('filtered PlanarOrbitData/',orbit_csv);
        [~,orbit_name,~] = fileparts(fullfile(folder_path,orbit_csv));
        output_file = fullfile(outputfolder,[orbit_name,'.mat']);
        fprintf(['Current family: ',orbit_name, '\n'])
        data = readmatrix(orbit_file);

        C = data(:,8);

        M = length(C);

        Wu_Section_Data = {};
        tu_Section_Data = {};
        Ws_Section_Data = {};
        ts_Section_Data = {};
        for j= 1:M
            fprintf(['Current family member: ',num2str(j), ' out of ', num2str(M) , '\n'])
            % Prepare to compute the unstable manifold
            %sgn = -1; % Set the sign for the manifold computation
            N = 50; % Numbeer of points on unstable manifold
            [Wup ,tup] = get_unst_manifold(orbit_file, C(j), N, 1,k);
            [Wum ,tum] = get_unst_manifold(orbit_file, C(j), N, -1,k);

            Wu_Section_Data{j} = [Wup;Wum];
            tu_Section_Data{j} = [tup;tum];

            % Computing corresponding stable manifold using time reversal symmetry
            if isempty(Wu_Section_Data{j})
                Ws_Section_Data{j} = [];
                ts_Section_Data{j} = [];
            else
                Wu = Wu_Section_Data{j};
                Ws = Wu;
                Ws(:,2:3) = -Wu(:,2:3);
                Ws_Section_Data{j} =  Ws;
                ts_Section_Data{j} = -tu_Section_Data{j};
            end
        end

    save(output_file, ...
        "Wu_Section_Data","tu_Section_Data","C" ...
        ,"Ws_Section_Data","ts_Section_Data",'-mat');
    end
end







%% --- Functions ---
function [X,PHI] = stm(tau,tspan,Yspan)
% Interpolates state transition matrix from propagation along single period
% of a PO.

T = tspan(end);
t = tau*T;
Y = interp1(tspan,Yspan(:,1:16),t);
X = interp1(tspan,Yspan(:,17:20),t)';
PHI = reshape(Y,[4,4]);

end

function [Wu, tu] = get_unst_manifold(orbit_file,C,N,sgn,section)
global mu
if isempty(gcp("nocreate"))
    parpool(14);
end
vareqn = @(t,x) var2D(t,x,mu);
varopt = odeset('RelTol',3e-10,'AbsTol',1e-10);

cr3bp = @(t,x) CR3BPMC2D(x,mu);

if section == 1
    cr3bp_opts = odeset('RelTol',3e-10,'AbsTol',1e-10, ...
        'Events',@(t,x) eX(t,x,mu));
elseif section ==2
    cr3bp_opts = odeset('RelTol',3e-10,'AbsTol',1e-10, ...
        'Events',@(t,x) mX(t,x,mu));

elseif section == 3
    cr3bp_opts = odeset('RelTol',3e-10,'AbsTol',1e-10, ...
        'Events',@(t,x) eY(t,x,mu));

elseif section == 4
    cr3bp_opts = odeset('RelTol',3e-10,'AbsTol',1e-10, ...
        'Events',@(t,x) mY(t,x,mu));
end
ep = 1e-5;
T_int = 20;

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

[~,ind] = max(real(diag(D)));
u = real(V(:,ind));


%Wsp_0 = zeros(4,N);
Wu_0 = zeros(4,N);



for i = 1 :length(tau)
    X = Yspan(i,17:20)';
    PHI = reshape(Yspan(i,1:16)',[4,4]);

    v = PHI*u;
    v = v/norm(v);

    if sgn == 1 && v(3) < 0
        v = -v;
    elseif sgn == -1 && v(3) > 0
        v = -v;
    end

    Wu_0(:,i) = X + ep*v;
end

Wu_temp = cell(1, N);
tu_temp = cell(1, N);

parfor i = 1:N
    [~, ~, te, We, ie] = ode45(cr3bp, [0,T_int],Wu_0(:,i),cr3bp_opts);
    Wu_temp{i} = We;
    tu_temp{i} = te;

end

Wu = vertcat(Wu_temp{:});
tu = vertcat(tu_temp{:});

end

%% Poincaré Sections
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
direction = 0;


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

direction = 0;

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


function [value, isterminal, direction] = x_axis_crossing(t,x,mu,pos)
% v = [x(3);x(4)];
% phi = acosd(v'*[1;0]);

r2 = sqrt((x(1)-1+mu)^2 + x(2)^2);
Rm = 1737/389703;


value(1) = x(2);
isterminal(1) = 0;
direction(1) = 0;

value(2) = r2 - Rm;
isterminal(2) = 1;
direction(2) = -1;
end

