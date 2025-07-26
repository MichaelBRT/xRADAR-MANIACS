clear, clc, close all

% global variables
global T_em D n TU_to_Days mu mass scalefactor
scalefactor = 1;

% plot colours
blue = [0.07,0.62,1.00];
orange = [0.988, 0.38, 0];
purple = [0.659, 0, 1];

% dimensional constants
T_em = 2.361 * 10^6;
D = 3.850*10^5;
n = 2*pi/T_em;

TU_to_Days = T_em/(2*pi*3600*24);

% Spacecraft mass
mass = 1000; % kg

% Discretization Points (Equal time spacing)
M = 1000;
tau = linspace(0,1,M);

% Order of Finite Fourier Series approximation
nffs  = 100;

% A matrix and derivatives
[A,Ad,Add] = getAmatrix2D(tau,nffs);


%% Computing departure and arrival orbits and their monodromy matrices
ind_i = 914;
ind_f = 933;


i_data = readmatrix("OrbitDataEarthMoon\Lyapunov (L1).csv");
f_data = readmatrix("OrbitDataEarthMoon\Lyapunov (L2).csv");

ind_i = find(i_data(:,1) == ind_i);
ind_f = find(f_data(:,1) == ind_f);

X0i = i_data(ind_i, [2,3,5,6])';
X0f = f_data(ind_f, [2,3,5,6])';

phi_0 = reshape(eye(4),[16,1]);

Y0i = [phi_0;X0i];
Y0f = [phi_0;X0f];

Ti = i_data(ind_i, 9);
Tf = f_data(ind_f, 9);

mu = 1.215058560962404E-2; 

vareqn = @(t,x) var2D(t,x,mu);
varopt = odeset('RelTol',3e-10,'AbsTol',1e-10);

[tspani,Yi] = ode113(vareqn,[0,Ti],Y0i,varopt);
[tspanf,Yf] = ode113(vareqn,[0,Tf],Y0f,varopt);

Mi = reshape(Yi(end,1:16),[4,4]);
Mf = reshape(Yf(end,1:16),[4,4]);

xi = Yi(:,17); yi = Yi(:,18);
xf = Yf(:,17); yf = Yf(:,18);


Re = 6378/389703;
Rm = 1737/389703;

xL1 = Lptpos(mu,1);
xL2 = Lptpos(mu,2);
xL3 = Lptpos(mu,3);


f1 = figure();
hold (gca, 'on')
plot(xi,yi,'LineWidth',2,'Color',blue)
plot(xf,yf,'LineWidth',2,'Color',purple);

draw_shaded_circle(gca,[mu,0],Re, [0.07,0.62,1.00])
draw_shaded_circle(gca,[1-mu,0],Rm, [0.1,0.1,0.1])
plot(xL1, 0, 'rx')
plot(xL2, 0, 'rx')
plot(xL3, 0, 'rx')
plot(1/2,sqrt(3)/2,'rx','LineWidth',1)
plot(1/2,-sqrt(3)/2,'rx','LineWidth',1)
hold (gca,'off')

grid on
axis equal

xlim([0.759 1.200])
ylim([-0.169 0.179])

%% Computing invariant manifolds
N = 200;
tau_span = linspace(0,1,N);

Wup_0 = zeros(4,N);
Wum_0 = zeros(4,N);

Wsp_0 = zeros(4,N);
Wsm_0 = zeros(4,N);


for i = 1:N
    [Wup_0(:,i),Wum_0(:,i)] = unst_IC(tau_span(i),tspani,Yi);
    [Wsp_0(:,i), Wsm_0(:,i)] = stab_IC(tau_span(i),tspanf,Yf); 
end

cr3bp = @(t,x) CR3BPMC2D(x,mu);
cr3bp_opts = odeset('RelTol',3e-10,'AbsTol',1e-10, ...
    'Events',@(t,x) Poincare1(t,x,mu));

T_int = 10*pi;

hold(gca,'on')
for i = 1:N
    [~,Wup_span,tup(i),Wup(i,:),~] = ode45(cr3bp, ...
        [0,T_int],Wup_0(:,i),cr3bp_opts);

    x = Wup_span(:,1);
    y = Wup_span(:,2);

    plot(x,y,'LineWidth',1,'Color',blue)

    [~,Wsm_span, tsm(i), Wsm(i,:),~] = ode45(cr3bp ...
        , [0,-T_int],Wsm_0(:,i),cr3bp_opts);
    
    x = Wsm_span(:,1);
    y = Wsm_span(:,2);

    plot(x,y,'LineWidth',1,'Color',purple)


end
hold(gca,'off')




diffs_Mat = zeros(N,N);

for i = 1:N
    for j = 1:N

        diffs_Mat(i,j) = sqrt(100*(Wup(i,2)-Wsm(j,2))^2 + ...
            (Wup(i,4) -  Wsm(j,4))^2);

    end
end

[minVal, linearIdx] = min(diffs_Mat(:));         
[rowIdx, colIdx] = ind2sub(size(diffs_Mat), linearIdx);  



%k = convhull(Wup(:,2),Wup(:,4));



%% Computing Poincare section of invariant manifolds
f2 = figure();
hold on
plot(Wup(:,2),Wup(:,4),'Color',blue,'LineStyle','none','Marker','.')
plot(Wup(rowIdx,2),Wup(rowIdx,4),'r*');
plot(Wsm(:,2),Wsm(:,4),'Color',purple,'LineStyle','none','Marker','.')
plot(Wsm(colIdx,2),Wsm(colIdx,4),'r*')
xlabel('y','FontSize',16,'Interpreter','latex')
ylabel('$\dot{y}$','FontSize',16,'Interpreter','latex')
hold off
grid on
% 
% g2 = figure();
% hold on
% plot(Wup(k,2),Wup(k,4))
% hold off
% grid on
% 



%%


cr3bp_opts = odeset('RelTol',3e-10,'AbsTol',1e-10);
XX0i = Wup(rowIdx,:)';
XX0f = Wsm(colIdx,:)';

% timespans of ballistic arcs
TTi = tup(rowIdx);
TTf = tsm(colIdx);

% span
[tti,XXi] = ode113(cr3bp,[0,-TTi],XX0i,cr3bp_opts);
[ttf,XXf] = ode113(cr3bp,[0,-TTf],XX0f,cr3bp_opts);


TOF = abs(TTi)+abs(TTf);

frac_span = linspace(0.01,0.2,200);

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



    %
    [ti,si] = ode45(cr3bp,[0,TOM],Si,cr3bp_opts);
    [tf,sf] = ode45(cr3bp,[0,-TOM],Sf,cr3bp_opts);

    tf = flip(tf) + TOM;
    sf = flip(sf);




    ssi = interp1(ti/TOM,si,tau);
    ssf = interp1(tf/TOM,sf,tau);


    S0 = interp_arc(ssi,ssf,M);

    x = S0(:,1);
    y = S0(:,2);

    xd = S0(:,1);
    yd = S0(:,2);


    [C,Cd, Cdd] = findCf2D(tau,Si,Sf, TOM);


    Cx = C(:,1);
    Cy = C(:,2);

    Cdx = Cd(:,1);
    Cdy = Cd(:,2);

    Cddx = Cdd(:,1);
    Cddy = Cdd(:,2);

    qx = pinv(A)*(x-Cx);
    qy = pinv(A)*(y-Cy);


    x = A*qx + Cx;
    y = A*qy + Cy;


    xd = 1/TOM*(Ad*qx + Cdx);
    yd = 1/TOM*(Ad*qy + Cdy);


    xdd = 1/TOM^2*(Add*qx + Cddx);
    ydd = 1/TOM^2*(Add*qy + Cddy);

    r1 = sqrt((x+mu).^2 + y.^2);
    r2 = sqrt((x-1+mu).^2 + y.^2);


    ux = xdd - 2*yd - x + (1-mu)*(x+mu)./r1.^3 + mu*(x-1+mu)./r2.^3;
    uy = ydd + 2*xd - y + (1-mu)*y./r1.^3 + mu*y./r2.^3;

    for i = 1:M
        U = [ux(i);uy(i)];
        a_nd(i) = U'*U;
    end



    % a = n^2*D*sqrt(a_nd);
    a = sqrt(a_nd);
    td = tau*TOM;

    Thr_max(j) = max(1000*a);
    dV(j) = n*D*trapz(td,a.^2);

end

[dVmin,dVmin_idx] = min(dV);
[Thrmin, Thrmin_idx] = min(Thr_max);



%%

t = tiledlayout(2,3);

nexttile();
plot(TOM_span*TU_to_Days,dV*1000,'LineWidth',2)
hold on
plot(TOM_span(Thrmin_idx)*TU_to_Days, dV(Thrmin_idx)*1000,'o','Color',orange,'MarkerFaceColor',orange)
plot(TOM_span(dVmin_idx)*TU_to_Days, dV(dVmin_idx)*1000,'o','Color',purple,'MarkerFaceColor',purple)
hold off
grid on
xlabel('Maneuver time $T_m$ [Days]','FontSize',16,'Interpreter','latex')
ylabel('Maneuver $\Delta V$ [m/s]','FontSize',16,'Interpreter','latex')
axis tight


nexttile(4);
plot(TOM_span*TU_to_Days,Thr_max*1000,'LineWidth',2)
hold on
plot(TOM_span(Thrmin_idx)*TU_to_Days, Thr_max(Thrmin_idx)*1000,'o','Color',orange,'MarkerFaceColor',orange)
plot(TOM_span(dVmin_idx)*TU_to_Days, Thr_max(dVmin_idx)*1000,'o','Color',purple,'MarkerFaceColor',purple)
hold off
grid on
xlabel('Maneuver time $T_m$ [Days]','FontSize',16,'Interpreter','latex')
ylabel('Max Thrust $\mathcal{T}$ [N]','FontSize',16,'Interpreter','latex')
axis tight

nexttile(2,[2,2]);
plot(Thr_max*1000,dV*1000,'LineWidth',2)
hold on
plot(Thr_max(Thrmin_idx)*1000, dV(Thrmin_idx)*1000,'o','Color',orange,'MarkerFaceColor',orange)
plot(Thr_max(dVmin_idx)*1000, dV(dVmin_idx)*1000,'o','Color',purple,'MarkerFaceColor',purple)
hold off
grid on
xlabel('Max Thrust $\mathcal{T}$ [N]','FontSize',16,'Interpreter','latex')
ylabel('Maneuver $\Delta V$ [m/s]','FontSize',16,'Interpreter','latex')
axis tight



% f3 = figure();
% hold (gca, 'on')
% plot(xi,yi,'LineWidth',2,'Color',blue)
% plot(xf,yf,'LineWidth',2,'Color',purple);
% 
% plot(XX0i(1),XX0i(2),'go','MarkerFaceColor','g')
% plot(XX0f(1),XX0f(2),'ro','MarkerFaceColor','r')
% 
% plot(XXi(:,1),XXi(:,2),'g','LineWidth',2)
% plot(XXf(:,1),XXf(:,2),'r','LineWidth',2)
% 
% plot(Si(1),Si(2),'*','Color',orange)
% plot(Sf(1),Sf(2),'*','Color',orange)
% 
% plot(si(:,1),si(:,2),'LineWidth',2,'Color',orange)
% plot(sf(:,1),sf(:,2),'LineWidth',2,'Color',orange)
% 
% plot(S0(:,1),S0(:,2),'k','LineWidth',2)
% 
% draw_shaded_circle(gca,[mu,0],Re, [0.07,0.62,1.00])
% draw_shaded_circle(gca,[1-mu,0],Rm, [0.1,0.1,0.1])
% plot(xL1, 0, 'rx')
% plot(xL2, 0, 'rx')
% plot(xL3, 0, 'rx')
% plot(1/2,sqrt(3)/2,'rx','LineWidth',1)
% plot(1/2,-sqrt(3)/2,'rx','LineWidth',1)
% hold (gca,'off')
% 
% 
% grid on
% axis equal
% 
% xlim([0.759 1.200])
% ylim([-0.169 0.179])
% 
% 
% 
% f4 = figure();
% 
% plot(tau,a*1000, 'LineWidth',2)
% grid on
% xlabel('$\tau$', 'FontSize',16,'Interpreter', 'latex')
% ylabel('$\mathcal{T} [kN]$', 'FontSize',16,'Interpreter', 'latex')

%% Optimization
mu = 1.215058560962404E-2;
TOM = TOM_span(Thrmin_idx);


% Generate initial guess for fmincon
% Initial and final states of maneuver arc
Si = interp1(tti,XXi,-TOM/2)';
Sf = interp1(ttf,XXf, TOM/2)';


[ti,si] = ode45(cr3bp,[0,TOM],Si,cr3bp_opts);
[tf,sf] = ode45(cr3bp,[0,-TOM],Sf,cr3bp_opts);

tf = flip(tf) + TOM;
sf = flip(sf);

ssi = interp1(ti/TOM,si,tau);
ssf = interp1(tf/TOM,sf,tau);


S0 = interp_arc(ssi,ssf,M);

x = S0(:,1);
y = S0(:,2);


[C,Cd, Cdd] = findCf2D(tau,Si,Sf, TOM);


Cx = C(:,1);
Cy = C(:,2);

Cdx = Cd(:,1);
Cdy = Cd(:,2);

Cddx = Cdd(:,1);
Cddy = Cdd(:,2);

qx = pinv(A)*(x-Cx);
qy = pinv(A)*(y-Cy);

Z0 = [qx;qy;log(scalefactor*TOM)];

% TOM0 = Z0(end);
% Tmin = TOM0 * 0.5;
% Tmax = TOM0 * 3.0;
% 
% lb = [-Inf(size(Z0,1)-1,1); Tmin];
% ub = [ Inf(size(Z0,1)-1,1); Tmax];

is_T_design = 1;

Tmax = 1;
params = struct('Wi',XXi,'Wf',XXf,'ti',tti,'tf',ttf, ...
    'M',M,'n',nffs,'A',A,'Ad',Ad,'Add',Add,'is_T_design',is_T_design ...
    ,'Thr_max',Tmax,'T',TOM);
%Z0 = [qx;qy];

%is_T_design = 0;

%Tmax = 1;
% trans_arc_params = struct('Wi',XXi,'Wf',XXf,'ti',tti,'tf',ttf, ...
%     'M',M,'n',nffs,'A',A,'Ad',Ad,'Add',Add,'is_T_design',is_T_design ...
%     ,'Thr_max',Tmax,'T',TOM);
% 
cost = @(Z) dVCostFun(Z,params);
cost2 = @(Z) max(thrust(Z,params));
cost3 = @(Z) combinedCost(Z,params);
constraint = @(Z) thrustConstraint(Z,params);

fmc_opts = optimoptions("fmincon","Display","iter-detailed" ...
    ,"EnableFeasibilityMode",true,"MaxFunctionEvaluations",5e3);
dummyfcn = @(Z) 0;

Z2 = fmincon(cost,Z0,[],[],[],[],[],[],constraint,fmc_opts);

% mmThr = Thrmin*1000;
% while Tmax-mmThr<0
%     Z1 = fmincon(cost2,Z0,[],[],[],[],[],[],constraint,fmc_opts);
%     Z2 = fmincon(cost,Z1,[],[],[],[],[],[],constraint,fmc_opts);
%     mmThr =  max(thrust(Z2,params));
%     Z0 = Z2;
% end
Z_opt = Z2;
qx = Z_opt(1:2*nffs+1);
qy = Z_opt(2*nffs+2:4*nffs+2);


% plotting stuff

[~,arc1]  = ode45(cr3bp,[0,TTi-TOM/2],XXi(end,:),cr3bp_opts);

% generating maneuver arc:
[C,Cd, Cdd] = findCf2D(tau,Si,Sf, TOM);

 Cx = C(:,1);    
 Cy = C(:,2);    

 x = A*qx + Cx;
 y = A*qy + Cy;

 arc2 = [x,y];

 [~,arc3]  = ode45(cr3bp,[0,TTf + TOM/2],XXf(end,:),cr3bp_opts);

 Thr = thrust(Z_opt,params);
 t_man = tau*TOM*TU_to_Days;


 f4 = figure();
 hold (gca, 'on')
 plot(arc1(:,1),arc1(:,2),'LineWidth',2,'Color',blue)
 plot(arc2(:,1),arc2(:,2),'LineWidth',2,'Color',orange)
 plot(arc3(:,1),arc3(:,2),'LineWidth',2,'Color',purple)

 plot(xi,yi,'LineWidth',2,'Color',blue)
 plot(xf,yf,'LineWidth',2,'Color',purple);

 draw_shaded_circle(gca,[mu,0],Re, [0.07,0.62,1.00])
 draw_shaded_circle(gca,[1-mu,0],Rm, [0.1,0.1,0.1])
 plot(xL1, 0, 'rx')
 plot(xL2, 0, 'rx')
 plot(xL3, 0, 'rx')
 plot(1/2,sqrt(3)/2,'rx','LineWidth',1)
 plot(1/2,-sqrt(3)/2,'rx','LineWidth',1)
 hold (gca,'off')

 grid on
 axis equal

 xlim([0.759 1.200])
 ylim([-0.169 0.179])


 %%
 
 f5 = figure();
 plot(t_man,Thr,'LineWidth',2,'Color',purple)
 grid on
 xlabel('maneuver time $t$ [days]','Interpreter','latex')
 ylabel('Thrust $\mathcal{T}$ [N]','Interpreter','latex')
 box on



%% functions
function J = combinedCost(Z,params)
    lambda = 0.75;
    dV = newcost(Z,params);
    T = max(thrust(Z,params));
    J = (1-lambda)*dV + lambda*T;        
end

function dV = dVCostFun(Z,params)
    dV = newcost(Z,params);
end

function T = thrust(Z,params)
global mass
    
    a = new_accel(Z,params);
    T = mass*a*1000; % N
end

function [c,ceq] = thrustConstraint(Z,params)
    T = thrust(Z,params);
    Tmax = params.Thr_max;
    c = T - Tmax;
    ceq = [];       
end

function [X,PHI] = stm(tau,tspan,Yspan)
% Interpolates state transition matrix from propagation along single period
% of a PO.

T = tspan(end);
t = tau*T;
Y = interp1(tspan,Yspan(:,1:16),t);
X = interp1(tspan,Yspan(:,17:20),t)';
PHI = reshape(Y,[4,4]);

end

function [xi_p, xi_m] = unst_IC(tau,tspan,Yspan)
    
    ep = 1e-6;
    T = tspan(end);
    X0 = Yspan(1,17:20)';

    [X,PHI] = stm(tau,tspan,Yspan);
    M = reshape(Yspan(end,1:16),[4,4]);

    [V,D] = eigs(M);

    [~,ind] = max(real(diag(D)));

    u = V(:,ind);
    u = PHI * u/norm(u);
    u = u/norm(u);

    xi_p = X + ep*u;
    xi_m = X - ep*u;
end

function [eta_p, eta_m] = stab_IC(tau,tspan,Yspan)
    ep = 1e-7;
    T = tspan(end);
    X0 = Yspan(1,17:20)';
    
    [X,PHI] = stm(tau,tspan,Yspan);
    M = reshape(Yspan(end,1:16),[4,4]);

    [V,D] = eigs(M);

    [~,ind] = min(real(diag(D)));
    u = V(:,ind);
    u = PHI*u/norm(u);
    u = u/norm(u);

    

    eta_p = X + ep*u;
    eta_m = X - ep*u;
end


function [value, isterminal, direction] = Poincare1(t,x,mu)
    v = [x(3);x(4)];
    phi = acosd(v'*[1;0]);

    if abs(phi) > 1
        value = x(1)-(1-mu);
    else
        value = [];
    end
    isterminal = 1;
    direction = 0;
end







%%
 function draw_shaded_circle(ax, center, radius, color)
        % draw_shaded_circle(center, radius, color)
        % Draws a filled circle at 'center' = [x, y] with given 'radius'
        % and optional 'color' (default is light blue).
        
        if nargin < 4
            color = [0.3, 0.6, 1]; % default color: light blue
        end
    
        theta = linspace(0, 2*pi, 100);  % 100 points around the circle
        u = center(1) + radius * cos(theta);
        v = center(2) + radius * sin(theta);
    
        fill(ax, u, v, color, 'EdgeColor', 'none');  % filled, no border
        %axis equal;  % keep aspect ratio correct
 end

function [A dA ddA] = getAmatrix(tau, n)
    % after merge
    m = length(tau);

    %%  A generation
    A = zeros(m, 1+2*n); 
    dA = zeros(m, 1+2*n);
    ddA = zeros(m, 1+2*n);


    for j = 1:m
        t = tau(j);

        % Row in A
        row = zeros(1, 2*n);
        row(1) = 0.5 * (1 - cos(2*pi*t));

        % Row in dA
        drow = zeros(1, 2*n);
        drow(1) = pi*sin(2*pi*t);

        % row in ddA
        ddrow = zeros(1, 2*n);
        ddrow(1) = 2*pi^2*cos(2*pi*t);


        for i = 3:(n + 2)
            idx = 2*(i - 3) + 2;
            ca = 0.5*((-1)^i - 1)*cos(pi*t) - ...
                0.5*((-1)^i + 1)*cos(2*pi*t) + cos(i*pi*t);

            dca = -0.5*pi*((-1)^i - 1)*sin(pi*t) +...
                pi*((-1)^i + 1)*sin(2*pi*t) - i*pi*sin(i*pi*t);
            
            ddca = -0.5*pi^2*((-1)^i - 1)*cos(pi*t) +...
                2*pi^2*((-1)^i+1)*cos(2*pi*t) - i^2*pi^2*cos(i*pi*t);



            cb = i/2*((-1)^i - 1)*sin(pi*t) - ...
                i/4*((-1)^i + 1)*sin(2*pi*t) + sin(i*pi*t);

            dcb = 0.5*i*pi*((-1)^i-1)*cos(pi*t) - ...
                0.5*i*pi*((-1)^i + 1) * cos(2*pi*t) +i*pi*cos(i*pi*t);

            ddcb = -0.5*i*pi^2*((-1)^i-1)*sin(pi*t) + ...
                i*pi^2*((-1)^i+1)*sin(2*pi*t) - i^2*pi^2*sin(i*pi*t);


            row(idx)   = ca;
            row(idx+1) = cb;

            drow(idx) = dca;
            drow(idx+1) = dcb;

            ddrow(idx) = ddca;
            ddrow(idx+1) = ddcb; 
        end
        A(j,:) = row;
        dA(j,:) = drow;
        ddA(j,:) = ddrow;
    end

end



 


 function dV = dVcostfcn(Z,tau,TOM,Si,Sf,A,Ad,Add,mu)
    T_em = 2.361 * 10^6;
    D = 3.850*10^5;
    n = 2*pi/T_em;

    M = length(tau);
    N = length(Z)/2;

    [C,Cd, Cdd] = findCf2D(tau,Si,Sf, TOM);


    Cx = C(:,1);
    Cy = C(:,2);

    Cdx = Cd(:,1);
    Cdy = Cd(:,2);

    Cddx = Cdd(:,1);
    Cddy = Cdd(:,2);

    qx = Z(1:N);
    qy = Z(N+1:end);


    x = A*qx + Cx;
    y = A*qy + Cy;


    xd = Ad*qx + Cdx;
    yd = Ad*qy + Cdy;


    xdd = Add*qx + Cddx;
    ydd = Add*qy + Cddy;

    r1 = sqrt((x+mu).^2 + y.^2);
    r2 = sqrt((x-1+mu).^2 + y.^2);


    ux = xdd - 2*yd - x + (1-mu)*(x+mu)./r1.^3 + mu*(x-1+mu)./r2.^3;
    uy = ydd + 2*xd - y + (1-mu)*y./r1.^3 + mu*y./r2.^3;

    for i = 1:M
        U = [ux(i);uy(i)];
        a_nd(i) = U'*U;
    end

    a = sqrt(a_nd);

    td = tau*TOM;

    % Thr_max = max(1000*a);
    dV = n*D*trapz(td,a.^2);
 end