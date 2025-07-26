clear, clc, close all
load("32_cycler_data.mat");
global mu

mu = 1.215058560962404E-2; 

filename = '(3,2) Cycler.csv';


% dimensional constants
T_em = 2.361 * 10^6;
D = 3.850*10^5;
n = 2*pi/T_em;

TU_to_Days = T_em/(2*pi*3600*24);

blue = [0.07,0.62,1.00];
orange = [0.988, 0.38, 0];
purple = [0.659, 0, 1];

Re = 6378/389703;
Rm = 1737/389703;

vareqn = @(t,x) var2D(t,x,mu);
varopt = odeset('RelTol',3e-10,'AbsTol',1e-10);


cr3bp = @(t,x) CR3BPMC2D(x,mu);

cr3bp_opts = odeset('RelTol',3e-10,'AbsTol',1e-10);
poincare_opts = odeset('RelTol',3e-10,'AbsTol',1e-10, ...
    'Events',@(t,x) Poincare1(t,x,mu));

x0 = dataset(1,1);
C0 = dataset(1,3);
T_days = dataset(1,2);
T = T_days/TU_to_Days;


C_span = linspace(3.1, C0, 1000);
C_span = flip(C_span);
orbit_family = zeros(length(C_span),11);
damping = 200;
converged = 0;
for i = 1:length(C_span)

    X0g = [x0,0,0,ydot_from_x(x0,C_span(i),-1)]';

    while ~converged
        damping = damping/2;
        [X0,T,M,converged] = SymDifCor_C(X0g,T/2-1,damping,0);
       % X0g = X0;
    end
    
    damping = 200;
    converged = 0;
    orbit_family(i,1) = 1;
    orbitfamily(i,2:3) = X0(1:2);
    orbitfamily(i,4) = 0;
    orbitfamily(i,5:6) = X0(3:4);
    orbitfamily(i,7) = 0; 
    orbitfamily(i,8) = C0; 
    orbitfamily(i,9) = T;
    orbitfamily(i,10) = stability_index(M);
    orbitfamily(i,11) = mu;

    x0 = X0(1);
end
familytable = array2table(orbitfamily,"VariableNames", ...
    ["Id","x0 (LU)", "y0 (LU)", "z0 (LU)", "vx0 (LU/TU)", "vy0 (LU/TU)" ...
    "vz0 (LU/TU)","Jacobi constant (LU2/TU2)", "Period (TU)", ...
    "Stability index", "Mass ratio"]);
writetable(familytable,filename);


% [t,X] = ode45(cr3bp,[0,T],X0,cr3bp_opts);
% 
% f2 = figure();
% plot(X(:,1),X(:,2),'LineWidth',2);
% grid on







function [value, isterminal, direction] = Poincare1(t,x,mu)
    tmin = 0.1;

    if t > tmin && x(1) > 1-mu
        value = x(2);
    else
        value = [];
    end
    isterminal = 0;
    direction = 0;

end


function ydot = ydot_from_x(x,C,pos)
global mu
% Computation of new ydot value
if nargin < 3
    pos = sign(C);
    if pos == -1
        C = -C ;
    end
end

% Gravitational parameter calculations
mu1=1-mu;
mu2=  mu;

U = -1/2*x.^2 - mu1./sqrt((x+mu2).^2)-mu2./sqrt((x-mu1).^2); 
ydot = (-C-2*U).^(1/2);

if pos==-1
    ydot =-ydot;
end

end

function I = stability_index(M)
    V  = eigs(M);
    [~,inds] = maxk(real(V)-1,2,'ComparisonMethod','abs');
    I = 1/2*sum(V(inds));
end