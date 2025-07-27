function [init_arc, U, thrust, deltaV_req,td,tTU, thrust_arc,XX0i,XX0f] = get_initial_arc(mu,app)





% Extract input parameters from interface
C1 = app.initJacobi.Value;
C2 = app.targetJacobi.Value;

orbit1_name = app.initDropdown.Value;
orbit2_name = app.targetDropdown.Value;

sectionID = app.sectionChoice.Value;



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

% Order of Finite Fourier Series approximation
nffs  = 100;

% A matrix and derivatives
[A,Ad,Add] = getAmatrix2D(tau,nffs);


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



if ~isempty(Wu)





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
    S0 = interp_arc(ssi,ssf,M);

    % Save the intermediate arc
    all_thrust_arcs{j} = S0;


    % thrust and DV computation

    % obtain x and y coordinates of intermediate arc
    x = S0(:,1);
    y = S0(:,2);

    % Compute C vectors (these contain boundary condition information)
    [C,Cd, Cdd] = findCf2D(tau,Si,Sf, TOM);


    Cx = C(:,1);
    Cy = C(:,2);

    Cdx = Cd(:,1);
    Cdy = Cd(:,2);

    Cddx = Cdd(:,1);
    Cddy = Cdd(:,2);

    % Compute vectors of Fourier coefficients
    qx = pinv(A)*(x-Cx);
    qy = pinv(A)*(y-Cy);


    % Compute approximate x and y coordinates from FS coeffs
    x = A*qx + Cx;
    y = A*qy + Cy;


    % Compute approximate xdot and ydot from FS coeffs
    xd = 1/TOM*(Ad*qx + Cdx);
    yd = 1/TOM*(Ad*qy + Cdy);


    % Compute approximate xddot and yddot from FS coeffs
    xdd = 1/TOM^2*(Add*qx + Cddx);
    ydd = 1/TOM^2*(Add*qy + Cddy);

    % From x,y and derivatives, compute thrust and dV using CR3BP dynamics
    r1 = sqrt((x+mu).^2 + y.^2);
    r2 = sqrt((x-1+mu).^2 + y.^2);


    ux{j} = xdd - 2*yd - x + (1-mu)*(x+mu)./r1.^3 + mu*(x-1+mu)./r2.^3;
    uy{j} = ydd + 2*xd - y + (1-mu)*y./r1.^3 + mu*y./r2.^3;

    for i = 1:M
        U = [ux{j}(i);uy{j}(i)];
        a_nd(i) = sqrt(U'*U);
    end



    % dimensionalize acceleration
    a = n^2*D*a_nd;   % km/s^2

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

thrust = [zeros(Nu,1);all_thrust{idx}';zeros(Ns,1)];

% Required delta V
deltaV_req = dV(idx);

% dimensional time corresponding to init_arc
td = (tau*TOF/n)'/(3600*24); % days
tTU = tau*TOF;               % TU   

else
    init_arc = [];
    tau = [];
    thrust = []; 
    deltaV_req = [];
    td = [];
    thrust_arc = [];
end



end



%% Functions

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


if abs(phi) > 10 && x(2) < -1737/389703 && abs(x(2))< 1
    value = x(1)-(1-mu);
else
    value = [];
end
isterminal = 1;
direction = 0;
end



function idx = findClosestJacobi(data, targetC)
    [~, idx] = min(abs(data(:,8) - targetC));
end


%% CR3BP ode function
function dx = CR3BPMC2D(x, mu)

    dx(1:2, 1) = x(3:4); %position derivatives = velocities
    
    
    r1 = sqrt(((x(1) + mu)^2) + (x(2)^2));
    r2 = sqrt(((x(1) - 1 + mu)^2) + (x(2)^2));
    
    dx(3, 1) = 2*x(4) + x(1) - (((1 - mu)*(x(1) + mu))/(r1^3)) - ((mu*(x(1) - 1 + mu))/(r2^3)); %xdtdt
    
    dx(4, 1) = -2*x(3) + x(2) - (((1 - mu)*x(2))/(r1^3)) - ((mu*x(2))/(r2^3)); %ydtdt

end


%% ~ Weighting Function  ~
%{       
Function to generate the weights used to get a rough arc between orbits.
Constraints:
1) w(0)=0 & w(1)=1
   -> a + b + c + d = 0.5

2) w'(0)=0 & w'(1)=0
   -> 7*a + 5*b + 3*c + d = 0

Program: xRADAR, Summer 2025
%}
% -------------------------------------------------------------------------
function [arc, w] = interp_arc( orbit1 , orbit2 , N) 


t = linspace(0, 1, N); % tau
% a = 0.02; % manually tuned
a = 0.0021416762037290; % optimally tuned
% b = 0.05; % manually tuned
b = 0.049040876142714; % optimally tuned
c = -2*b - 3.25*a - 0.25;
d = b + 2.25*a + 0.75;

tau = 2*t - 1;
w = a * tau.^7 + b * tau.^5 + c * tau.^3 + d * tau + 0.5;
w = min(w, 1);  % Soft clip anything above 1
w = w';  % ensures w is N x 1
arc = (1 - w) .* orbit1 + w .* orbit2;
% arc = (1-w)*orbit1 + w*orbit2;
  

end



function [A, dA, ddA] = getAmatrix(tau, n)

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


function [Cf,Cfd,Cfdd] = findCf2D(tou, Xi, Xf, TOM)
    m = length(tou);
    Cf = zeros(m,2);
    Cfd = zeros(m,2);
    Cfdd = zeros(m,2);

    
    for index = 1:2
        % Boundary conditions
        f0 = Xi(index); f1 = Xf(index);
        df0 = TOM*Xi(index+2); df1 = TOM*Xf(index+2);
        

        for j = 1:m
            t = tou(j);
            Cf(j,index) = 0.5 * (f0 - f1) * cos(pi*t) + ...
                0.5/pi * (df0 - df1) * sin(pi*t) + ...
                0.5 * (f0 + f1) * cos(2*pi*t) + ...
                0.25/pi * (df0 + df1) * sin(2*pi*t);
            Cfd(j,index) = -(pi/2) * (f0-f1) * sin(pi*t) + ...
            (df0-df1) * cos(pi*t) * 0.5 -pi*(f0+f1) * sin(2*pi*t) + ...
            (df0+df1)*cos(2*pi*t) * 0.5;

            Cfdd(j,index) = -(pi*pi*0.5)*(f0-f1)*cos(pi*t) - ...
            0.5 * pi * (df0-df1) * sin(pi*t) - ...
            2*  pi * pi * (f0+f1) * cos(2*pi*t) - ...
            pi * (df0+df1) * sin(2*pi*t);

        end
    end
end