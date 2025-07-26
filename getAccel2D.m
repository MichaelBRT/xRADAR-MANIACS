function a =  getAccel2D(Z,params)
global D n mu scalefactor

check = params.is_T_design;
if check
    T = exp(Z(end))/scalefactor;
else
    T = params.T;
end

M = params.M;
nffs = params.n;

N = 2*nffs+1;

tau = linspace(0,1,M);


Wi = params.Wi;
Wf = params.Wf;

ti = params.ti; 
tf = params.tf; 


Si = interp1(ti, Wi, -T/2);
Sf = interp1(tf, Wf, T/2);


qx = Z(1:N);
qy = Z(N+1:2*N);

A = params.A; Ad = params.Ad; Add = params.Add;
[C,Cd, Cdd] = findCf2D(tau,Si,Sf,T);


Cx = C(:,1);    Cdx = Cd(:,1);    Cddx = Cdd(:,1);
Cy = C(:,2);    Cdy = Cd(:,2);    Cddy = Cdd(:,2);

x = A*qx + Cx;  xd = 1/T * (Ad*qx + Cdx); xdd = 1/T^2 * (Add*qx + Cddx);
y = A*qy + Cy;  yd = 1/T * (Ad*qy + Cdy); ydd = 1/T^2 * (Add*qy + Cddy);

r1 = sqrt((x+mu).^2 + y.^2);
r2 = sqrt((x-1+mu).^2 + y.^2);

ux = xdd - 2*yd - x + (1-mu)*(x+mu)./r1.^3 + mu*(x-1+mu)./r2.^3;
uy = ydd + 2*xd - y + (1-mu)*y./r1.^3 + mu*y./r2.^3;


for i = 1:M
    U = [ux(i);uy(i)];
    a_nd(i) = sqrt(U'*U);
end

a = n^2*D*a_nd;



end

