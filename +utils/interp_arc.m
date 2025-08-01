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
function [arc, Darc] = interp_arc(orbit1 , orbit2 , N, T, mu) 
cr3bp = @(X) utils.pcr3bp(0,X,mu);

for i = 1 : N
    dorbit1(i,:) = cr3bp(orbit1(i,:)')';
    dorbit2(i,:) = cr3bp(orbit2(i,:)')';
end

t = linspace(0, 1, N); % tau
% a = 0.02; % manually tuned
a = -0.583333333333333; % optimally tuned
% b = 0.05; % manually tuned
b = 0.9166666666666667; % optimally tuned
c = -3*a - 2*b - 0.25;
d = 2*a + b + 0.75;

tau = 2*t - 1;
w = a * tau.^7 + b * tau.^5 + c * tau.^3 + d * tau + 0.5;
dw = 2*(7*a*tau.^6 + 5*b*tau.^4 + 3*c*tau.^2 + d);
w = min(w, 1);  % Soft clip anything above 1
w = w';  % ensures w is N x 1
dw = dw';
arc = (1 - w) .* orbit1 + w .* orbit2;
Darc = 1/T*dw.*(orbit2 - orbit1) + (1-w).*dorbit1 + w.*dorbit2;

% arc = (1-w)*orbit1 + w*orbit2;
  


%{
% Anonymous function to optimize
opt_fun = @(x) max(x(1)*tau.^7 + x(2)*tau.^5 + ...
                  (-2*x(2) - 3.25*x(1) - 0.25)*tau.^3 + ...
                  (x(2) + 2.25*x(1) + 0.75)*tau + 0.5) - 1;

x0 = [0.02, 0.05];  % initial guess for [a, b]
x_opt = fsolve(opt_fun, x0)
%}

%{
% Plot
plot(t, w, 'LineWidth', 2);
xlabel('\tau')
ylabel('w(\tau)')
title(sprintf('Weighting Function w(\\tau) for a = %.2f', a))
ylim([0 1]);
xlim([0 1]);
grid on
%}


end
