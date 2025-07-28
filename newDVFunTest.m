clear, clc, close all


mu = 1.215058560962404E-2;
m = 1000;
blue = [0.07,0.62,1.00];
orange = [0.988, 0.38, 0];
purple = [0.659, 0, 1];


% ode function
tol = 1e-10;
cr3bp = @(t,x) utils.pcr3bp(t,x,mu);


% Initial parameters (same as in main, this can be trimmed for efficiency)
T_em = 2.361 * 10^6;
D = 3.850*10^5;
n = 2*pi/T_em;
TU_to_Days = T_em/(2*pi*3600*24);


M = 1000;
tau = linspace(0,1,M);


mass = 1000;

C1 = 3.17;
C2 = 3.15;


sectionID = 'mY';
orbit1_name = 'Lyapunov (L1).csv';
orbit2_name = 'Lyapunov (L2).csv';

orbit1_file = fullfile('filtered PlanarOrbitData/',orbit1_name);
orbit2_file = fullfile('filtered PlanarOrbitData/',orbit2_name);

i_data = readmatrix(orbit1_file);
f_data = readmatrix(orbit2_file);

ind_i = utils.findClosestJacobi(i_data,C1);
ind_f = utils.findClosestJacobi(f_data,C2);

[~,orbit1_name,~] = fileparts(orbit1_file);
[~,orbit2_name,~] = fileparts(orbit2_file);

orbit1_Poincare_File = fullfile('Poincaré Section Data/', ...
    [sectionID,'/',orbit1_name, '.mat']);
orbit2_Poincare_File = fullfile('Poincaré Section Data/', ...
    [sectionID,'/',orbit2_name, '.mat']);

% Load Poincaré section data for both orbits
load(orbit1_Poincare_File);

Wu = Wu_Section_Data{ind_i};
tu = tu_Section_Data{ind_i};

load(orbit2_Poincare_File);
Ws = Ws_Section_Data{ind_f};
ts = ts_Section_Data{ind_f};



% Search for "closest" crossing
diffs_Mat = zeros(size(Wu,1),size(Ws,1));



for i = 1:size(Wu,1)
    for j = 1:size(Ws,1)

        diffs_Mat(i,j) = sqrt((Wu(i,2)-Ws(j,2))^2 +...
            0.1*(Wu(i,4)-Ws(j,4))^2 + 0.1*(Wu(i,3)-Ws(j,3))^2);

    end
end
%
[~, linearIdx] = min(diffs_Mat(:));      
[rowIdx, colIdx] = ind2sub(size(diffs_Mat), linearIdx);  

% rowIdx - index of crossing on the unstable manifold
% colIdx - index of crossing on the stable manifold 



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

TOM = 0.09*TOF;

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


% Extract initials state of orbits in question
X0i = i_data(ind_i, [2,3,5,6])';
X0f = f_data(ind_f, [2,3,5,6])';

%
Ti = i_data(ind_i, 9);
Tf = f_data(ind_f, 9);


cr3bp_opts = odeset('RelTol',3e-13,'AbsTol',1e-13);

[~,Xi] = ode45(cr3bp,[0,Ti],X0i,cr3bp_opts);
xi = Xi(:,1); yi = Xi(:,2);
[~,Xf] = ode45(cr3bp,[0,Tf],X0f,cr3bp_opts);
xf = Xf(:,1); yf = Xf(:,2);

f1 = figure();
f1_ax = gca;
hold (gca, 'on')


plot(xi,yi,'LineWidth',2,'Color',blue)
plot(xf,yf,'LineWidth',2,'Color',purple);

plot(ssi(:,1),ssi(:,2),'LineWidth',1.5,'Color','r')
plot(ssf(:,1),ssf(:,2),'LineWidth',1.5,'Color','g')

plot(XXi(:,1),XXi(:,2),'LineWidth',1.5,'Color','r')
plot(XXf(:,1),XXf(:,2),'LineWidth',1.5,'Color','g')


utils.drawEarthMoonSystem(f1_ax,0,C1)
xlim([0.7,1.3]);
ylim([-0.3,0.3])


a = 0.0021416762037290;
b = 0.049040876142714;

[~,dV,thrust,~,c] = costfun_WO5(a,b,TOM,ssi,ssf,mu,m);
dV = dV*D*n*1000

hold(f1_ax,"on")
plot(f1_ax,c(:,1),c(:,2),'LineWidth',2,'Color',orange)

utils.theoretical_min_dV(i_data(ind_i,8),f_data(ind_f,8),mu)

num_pts = 25;



a = linspace(-1,1,num_pts);
b = linspace(-1,1,num_pts);

[aa,bb] = meshgrid(a,b);

all_curves = {};

num_TOM = 100;
TOM = linspace(0.05,0.5,num_TOM)*TOF;
parfor k = 1:num_TOM
    dV = zeros(num_pts,num_pts);
    for i = 1:num_pts
        for j=  1:num_pts
            a = aa(i,j); b = bb(i,j);
            [~,dV(i,j),~,~,~] = costfun_WO5(a,b,TOM(k),ssi,ssf,mu,m);


        end
    end
    [minVal, linearIdx] = min(dV(:));           % Linear index of min
    [row, col] = ind2sub(size(dV), linearIdx);  % Convert to (row, col)
    dVk(k) = minVal*D*n*1000;
    ak(k) = aa(row,col);
    bk(k) = bb(row,col);
    ck(k) = -3*ak(k) - 2*bk(k) -0.25;
    dk(k) = 2*ak(k) + bk(k) + 0.75;

end




f2 = figure();
plot(TOM,ck,'LineWidth',2);
grid on
% contourf(aa, bb, dV*D*n*1000, 50, 'LineColor', 'none');
% axis equal;
% colorbar;
% title('Contourf heatmap of \Delta V(a,b)');
% xlabel('a'); ylabel('b');
