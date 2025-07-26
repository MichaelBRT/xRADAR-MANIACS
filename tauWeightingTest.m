clear, clc, close all;

addpath("orbit database\")

mu = 0.01215;
orbitdata = readmatrix("Halo (L2) (Southern).csv");
[~,ind] = min(abs(orbitdata(:,8)-3.07672711));

T = orbitdata(ind,9);
X0 = orbitdata(ind,2:7)';

opts = odeset('RelTol',1e-10,'AbsTol',1e-11);
odefcn = @CR3BPMC;

[t,X] =  ode45(odefcn,[0,T],X0,opts);

x = X(:,1); y = X(:,2) ; z = X(:,3);
Rm = 1737/389703;
[Xm,Ym,Zm] = ellipsoid(1,0,0,Rm,Rm,Rm,400);



% interpolating points:

for i = 1:length(t)
    a(i) = grav_accel(X(i,:),mu);
end

I = trapz(t/T,a);
npts = 30;

Fvals = cumtrapz(t/T, a);

target_vals = linspace(0,1,npts+1);

tau = interp1(Fvals,t/T,target_vals);

xunif = interp1(t/T,X(:,1),tau);
yunif = interp1(t/T,X(:,2),tau);
zunif = interp1(t/T,X(:,3),tau);



f = figure('Position',[488   242   210   420]);
plot3(x,y,z,'LineWidth',2)
hold on
plot3(xunif,yunif,zunif,'ro','MarkerFaceColor','r','LineStyle','none')
surf(Xm,Ym,Zm,'EdgeColor',[0.1,0.1,0.1],'FaceColor','flat')
grid on
axis equal

f2 = figure();
plot(t/T, a, 'LineWidth',2)
grid on


function a = grav_accel(X,mu)
    x = X(1); y = X(2); z = X(3);

    r1 = [x+mu,y,z];
    r2 = [x-(1-mu),y,z];

    a = (1-mu)*r1/norm(r1)^3 + mu*r2/norm(r2)^3;
    a = norm(a);
end
