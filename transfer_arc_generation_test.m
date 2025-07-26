clear, clc, close all;

tic;

T_em = 2.361 * 10^6;
D = 3.850*10^5;
n = 2*pi/T_em;
mass = 1000; %kg

nd_to_N = D*n^2*mass*1000;

TU_to_Days = T_em/(2*pi*3600*24);

% Spacecraft mass
mass = 1000; % kg

C_i = 3.17;
C_f = 3.16;
orbit1_file = fullfile("filtered PlanarOrbitData/Lyapunov (L1).csv");
orbit2_file = fullfile("filtered PlanarOrbitData/Lyapunov (L2).csv");

filename_init = fullfile("Poincaré Section Data/mY/Lyapunov (L1).mat");
filename_final = fullfile("Poincaré Section Data/mY/Lyapunov (L2).mat");

i_data = readmatrix(orbit1_file);
f_data = readmatrix(orbit2_file);

idx1 = utils.findClosestJacobi(i_data,C_i);
idx2 = utils.findClosestJacobi(f_data,C_f);

mu = f_data(idx1,end);



load(filename_init);
load(filename_final);

Wu_sec = Wu_Section_Data{idx1};
Ws_sec = Ws_Section_Data{idx2};

tu_sec = tu_Section_Data{idx1};
ts_sec = ts_Section_Data{idx2};

[rowIdx,colIdx] = utils.get_manifold_intersection(Wu_sec,Ws_sec);

f1 = figure(2);
f1.Theme = 'dark';
hold on

scatter(Wu_sec(:,2),Wu_sec(:,4),7.5,'r','filled')
plot(Wu_sec(rowIdx,2),Wu_sec(rowIdx,4),'Color','r','Marker','x','LineWidth',1.5 ...
    ,'MarkerSize',10)
plot(Ws_sec(colIdx,2),Ws_sec(colIdx,4),'Color','g','Marker','x','LineWidth',1.5 ...
    ,'MarkerSize',10)
scatter(Ws_sec(:,2),Ws_sec(:,4),7.5,'g','filled')

grid on

Xu = Wu_sec(rowIdx,:)';
Xs = Ws_sec(colIdx,:)';

tu = tu_sec(rowIdx);
ts = ts_sec(colIdx);

tol = 1e-13;
cr3bp = @(t,x) utils.pcr3bp(t,x,mu);
cr3bp_opt = odeset('RelTol',3*tol,'AbsTol',tol);

[Tu,Wu] = ode113(cr3bp,[2*tu,0],Xu,cr3bp_opt);
Tu = flip(Tu); Wu = flip(Wu);
[Ts,Ws] = ode113(cr3bp,[2*tu,2*tu-1.5*ts],Xs,cr3bp_opt);
TOF  =2*tu - 1.5*ts;

m = 1000; % number of DPs
nffs = 100;  % order of FFS

nn = 50; % number of thrust arcs to try
frac = linspace(0.05,0.1,nn);

tol = 1e-10;
cr3bp_opt = odeset('RelTol',3*tol,'AbsTol',tol);

AAA = utils.getAmatrix(m,nffs);

for i = 1:nn
    TOM(i) = frac(i)*TOF;
    X0 = (interp1(Tu,Wu,2*tu-TOM(i)/2,"spline"))';
    [~,Xminus] = ode113(cr3bp,linspace(0,TOM(i),m),X0,cr3bp_opt);
    
    Xf = (interp1(Ts,Ws,2*tu + TOM(i)/2,"spline"))';
    [~,Xplus] = ode113(cr3bp,linspace(0,-TOM(i),m),Xf,cr3bp_opt);
    Xplus = flip(Xplus);

    X{i} = utils.interp_arc(Xminus,Xplus,m);
    CCC{i} = utils.getCf(linspace(0,1,m),X0,Xf,TOM(i));

    [~,~,u{i}] = getControls_from_State(X{i},TOM(i),AAA,CCC{i},mu);
    umax(i) = max(u{i});
end

[thrmax_nd,minThr_idx] = min(umax);

thrmax = thrmax_nd*nd_to_N;
dV = trapz(linspace(0,1,m)*TOM(minThr_idx)/n,u{minThr_idx}*D*n^2*1000);

annotation_text = {['$\Delta V = $', num2str(dV) 'm/s'],...
    ['$\max \mathcal{T} = $' num2str(thrmax), 'N']};

f2 = figure();
ax = gca;
f2.Theme = 'dark';

%
ax.Units = 'normalized';
ax_pos = ax.Position;

% Compute top-left of axes
x_left = ax_pos(1);
y_top  = ax_pos(2) + ax_pos(4);


hold(ax,"on")
plot(Wu(:,1),Wu(:,2),'r','LineWidth',2)
plot(Ws(:,1),Ws(:,2),'g','LineWidth',2)
plot(X{minThr_idx}(:,1),X{minThr_idx}(:,2),'Color','y','LineWidth',2)
utils.drawEarthMoonSystem(ax,1,C_i)
xlim([0.8,1.2])
ylim([-0.2,0.2])
%
ax.Units = 'normalized';

% Get full figure-space position of plot box (not just axes container)
xl = xlim;
yl = ylim;
text(xl(1), yl(2), annotation_text, ...
     'HorizontalAlignment', 'left', ...
     'VerticalAlignment', 'top', ...
     'FontSize', 12, ...
     'Interpreter','latex', ...
     'Color','y')

toc;

function [m,n] = allAMatrices2Dim(allAMatrices)
    N = size(allAMatrices,2);
    n = 1/2*(N-1);

    m = 1/3*size(allAMatrices,1);
end

function [ux,uy,u] = getControls_from_State(X,TOM,AAA,CCC,mu)
    
    [m,~] = allAMatrices2Dim(AAA);
    A = AAA(1:m,:);
    dA = AAA(m+1:2*m,:);
    ddA = AAA(2*m+1:3*m,:);

    % Extract boundary information
    Cx = CCC(:,1); Cy = CCC(:,2);
    dCx = CCC(:,3); dCy = CCC(:,4);
    ddCx = CCC(:,5); ddCy = CCC(:,6);

    % Compute FFS coefficients for the given thrust arc:
    qx = pinv(A)*(X(:,1)-Cx);
    qy = pinv(A)*(X(:,2)-Cy);

    % x and y cords
    x = A*qx + Cx;
    y = A*qy + Cy;

    % x and y velocities
    vx = 1/TOM*(dA*qx + dCx);
    vy = 1/TOM*(dA*qy + dCy);

    % x and y accelerations
    ax = 1/TOM^2*(ddA*qx + ddCx);
    ay = 1/TOM^2*(ddA*qy + ddCy);

   
    % distance from primaries
    r1 = sqrt((x+mu).^2 + y.^2);
    r2 = sqrt((x-1+mu).^2 + y.^2);

    % Control inputs in x and y
    ux = ax - 2*vy - x + (1-mu)*(x+mu)./r1.^3 + mu*(x-1+mu)./r2.^3;
    uy = ay + 2*vx - y + (1-mu)*y./r1.^3 + mu*y./r2.^3;

    % Control input magnitude
    u = sqrt(ux.^2 + uy.^2);

end

function u =  getControls_from_Design(Z,TOM,AAA,CCC,mu)
    N = length(Z);

    
    qx = Z(1:N/2);
    qy = Z(N/2+1:N);

    [m,~] = allAMatrices2Dim(AAA);

    A = AAA(1:m,:);
    dA = AAA(m+1:2*m,:);
    ddA = AAA(2*m+1:3*m,:);
    
    Cx = CCC(:,1); Cy = CCC(:,2);
    dCx = CCC(:,3); dCy = CCC(:,4);
    ddCx = CCC(:,5); ddCy = CCC(:,6);

    % x and y cords
    x = A*qx + Cx;
    y = A*qy + Cy;

    % x and y velocities
    vx = 1/TOM*(dA*qx + dCx);
    vy = 1/TOM*(dA*qy + dCy);

    % x and y accelerations
    ax = 1/TOM^2*(ddA*qx + ddCx);
    ay = 1/TOM^2*(ddA*qy + ddCy);

    % distance from primaries
    r1 = sqrt((x+mu).^2 + y.^2);
    r2 = sqrt((x-1+mu).^2 + y.^2);

    % Control inputs in x and y
    ux = ax - 2*vy - x + (1-mu)*(x+mu)./r1.^3 + mu*(x-1+mu)./r2.^3;
    uy = ay + 2*vx - y + (1-mu)*y./r1.^3 + mu*y./r2.^3;

    % Control input magnitude
    u = sqrt(ux.^2 + uy.^2);
end


function dV =  getDV_from_Design(Z,TOM,AAA,CCC,mu)

    T_em = 2.361 * 10^6;
    D = 3.850*10^5;
    n = 2*pi/T_em;
    mass = 1000; %kg

    [m,~] = allAMatrices2Dim(AAA);
    nd_to_N = D*n^2*mass*1000;

    thrust = getControls_from_Design(Z,TOM,AAA,CCC,mu)*nd_to_N;
    time = linspace(0,1,m)*TOM/n; % s

    dV = 1/mass*trapz(time,thrust); % m/s
end

function [c,ceq] = low_thrust_constraint(Z,TOM,AAA,CCC,mu,thrustlim)
    T_em = 2.361 * 10^6;
    D = 3.850*10^5;
    n = 2*pi/T_em;
    mass = 1000; %kg

    nd_to_N = D*n^2*mass*1000;
    thrust = getControls_from_Design(Z,TOM,AAA,CCC,mu)*nd_to_N;

    c = thrust.^2 - thrustlim.^2;
    ceq = [];
end


function stop = designOutput(Z,~,state,TOM,AAA,CCC,mu,ax)
    stop = 0;
    switch state
        case 'iter'
              cla(ax);
              drawnow;
              utils.drawEarthMoonSystem(ax,1);
              drawTransferFromDesign(Z,TOM,AAA,CCC,mu,ax);
              drawnow;
              % Make updates to plot or guis as needed
        case 'init'

              % Setup for plots or guis
              cla(ax);
              drawnow;
              utils.drawEarthMoonSystem(ax,1);
              drawTransferFromDesign(Z,TOM,AAA,CCC,mu,ax);
              drawnow;

        case 'done'
              % Cleanup of plots, guis, or final plot
              cla(ax);
              drawnow;
              utils.drawEarthMoonSystem(ax,1);
              drawTransferFromDesign(Z,TOM,AAA,CCC,mu,ax);
              drawnow;
    otherwise
    end
end

function drawTransferFromDesign(Z,TOM,AAA,CCC,mu,ax)
    orange = [0.988, 0.38, 0];

    T_em = 2.361 * 10^6;
    D = 3.850*10^5;
    n = 2*pi/T_em;
    mass = 1000; %kg

    nd_to_N = D*n^2*mass*1000;
    thrust = getControls_from_Design(Z,TOM,AAA,CCC,mu)*nd_to_N;
    thrmax = max(thrust);

    N = length(Z);

    qx = Z(1:N/2);
    qy = Z(N/2+1:N);

    [m,~] = allAMatrices2Dim(AAA);

    A = AAA(1:m,:);    
    Cx = CCC(:,1); Cy = CCC(:,2);


    % x and y cords
    x = A*qx + Cx;
    y = A*qy + Cy;

    dV = getDV_from_Design(Z,TOM,AAA,CCC,mu);
    
    
    hold(ax,"on");
    plot(ax,x,y,'LineWidth',2,"Color",orange);
    plot(ax,x(1),y(1),'o','LineWidth',1.5,'Color',orange, ...
        'MarkerFaceColor','w');
    plot(ax,x(end),y(end),'o','LineWidth',1.5,'Color',orange, ...
        'MarkerFaceColor','w');
    hold(ax,"on");
    %
    ax.Units = 'normalized';
    ax_pos = ax.Position;

    % Compute top-left of axes
    x_left = ax_pos(1);
    y_top  = ax_pos(2) + ax_pos(4);

    annotation_text = {['$\Delta V = $', num2str(dV) 'm/s'],...
        ['$\max \mathcal{T} = $' num2str(thrmax), 'N']};
    
    % Place annotation textbox at top-left corner of axes
    annotation('textbox', [x_left, y_top - 0.05, 0.1, 0.05], ...
        'String', annotation_text, ...
        'FitBoxToText', 'on', ...
        'VerticalAlignment', 'top', ...
        'EdgeColor', 'none', ...
        'Interpreter','latex');

end

