function [X0,T, M, converged] = SymDifCor_C(X0g,tmin, damping, showplot)
    global mu
    vareqn = @(t,x) var2D(t,x,mu);
    varopt = odeset('RelTol',3e-10,'AbsTol',1e-10, ...
        'Events',@(t,x) xaxiscrossing(t,x,tmin),'OutputFcn',[]);
    
    
    cr3bp = @(t,x) CR3BPMC2D(x,mu);
    % cr3bp_opts = odeset('RelTol',3e-12,'AbsTol',1e-12, ...
    %     'Events',@(t,x) Poincare1(t,x,mu));

    iter_max = 1000;
    tol = 1e-10;
    count = 0;
    vx = 1;


    C = Jconst(X0g);
    pos = sign(X0g(4));
            
    x0 = X0g(1);
    vy0 = ydot_from_x(x0,C,pos);
    PHI0 = reshape(eye(4),[16,1]);
    % if showplot
    %     f = figure();
    %     grid on
    % end

    while count < iter_max
        count = count+1;

        [t,PHI_vec] = ode45(vareqn,[0,10*pi],[PHI0;X0g],varopt);
        % if showplot
        %     X = PHI_vec(:,17);
        %     Y = PHI_vec(:,18);
        % 
        %     plot(X,Y,'LineWidth',2)
        %     grid on
        %     drawnow;
        % end

        vx = PHI_vec(end,19);
        if mod(count,10) == 0
            fprintf('Attempt %3d: vx = %.3e\n', count, vx);
            %pause(0);
        end
        if abs(vx) < tol
            fprintf('Attempt %3d converged. Exiting Loop...\n',count)
            converged = 1;
            break;
        end

        x = PHI_vec(end,17);
        y = PHI_vec(end,18);

        vy = PHI_vec(end,20);

        P = reshape(PHI_vec(end,1:16),[4,4]);
        
        Ux0 = -x0 + (1-mu)*(x0+mu)/abs(x0+mu)^3 +...
        mu*(x0-1+mu)/abs(x0-1+mu)^3;
        
        
        r1 = sqrt((x+mu)^2+y^2);
        r2 = sqrt((x-1+mu)^2+y^2);

        Ux1 = -x + (1-mu)*(x+mu)/r1^3 + mu*(x-1+mu)/r2^3;
        ax = 2*vy - Ux1;


        dx0 = vx*vy/(ax*P(2,1))/(1-P(2,4)/P(2,1)*1/vy0*Ux0- ...
            vy/ax*1/P(2,1)*(P(3,1)-P(3,4)*1/vy0*Ux0));

        x0 = x0 + damping*dx0;

        vy0 = ydot_from_x(x0,C,pos);

        X0g = [x0,0,0,vy0]';

        
        

    end

    if count >= iter_max
        converged = 0;
        fprintf('Maximum number of iterations reached.')
    end

    X0 = X0g;
    T = 2*t(end);
    G = diag([1,-1,-1,1]);
    M = G * (P \ G) * P;

end


function PHIdot = var2D(~,PHI,mu)
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

function [value,isterminal,dir] = xaxiscrossing(t,x,tmin)

if t > tmin
    value = x(18);
else
    value = [];
end
isterminal = 1;
dir = 0;
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

function C = Jconst(X)
    % Function calculates Jacobi constant given a state
    global mu

    x = X(1); 
    y = X(2);
    xdot = X(3); 
    ydot = X(4);

   
    mu1=1-mu;
    mu2=  mu;

    r1 = ((x+mu2)^2 + y^2)^(1/2);
    r2 = ((x-mu1)^2 + y^2)^(1/2);
    
    U = -1/2*(x^2+y^2)-mu1/r1-mu2/r2; %-1/2*mu1*mu2;
    
    C = -(xdot^2 + ydot^2)-2*U;
end