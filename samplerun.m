

m = 60;
n_terms = 5;                              
tau = linspace(0, 1, m)'; 

Xx = 0.5 * rand(2*n_terms+1, 3);  
Xi = [0; 0; 0; 0; 0; 0];             
Xf = [1; 1; 1; 0; 0; 0];            

[Cf, Cfd, Cfdd] = findCf(tau, Xi, Xf,m);
[A, Ad, Add] = getAmatrix(tau, n_terms);

x     = A * Xx + Cf;
xdot  = Ad * Xx + Cfd;
xddot = Add * Xx + Cfdd;
figure;
subplot(3,1,1)
plot(tau, x, 'LineWidth', 1.5); grid on;
title('x');

subplot(3,1,2)
plot(tau, xdot, 'LineWidth', 1.5); grid on;
 title('xd');

subplot(3,1,3)
plot(tau, xddot, 'LineWidth', 1.5); grid on;
title('xdd');
