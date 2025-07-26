% function [cf, ax, Xx, x] = findmatrices(f, tou, acoeff, bcoeff)
% nffs = length(acoeff);
% N = length(tou);

function [A, Cf] = buildmatrices(tau, Xi, Xf, nffs)
% Construct the Ax matrix, Cf constant vector, and evaluate x and X
% tau: (m x 1) vector of normalised time in [0, 1]
% fi, ff: initial and final position
% dfi, dff: initial and final velocity
% acoeff, bcoeff: Fourier series coefficients
% index:
%   1 - x
%   2 - y
%   3 - z
% hello JMS 
    m = length(tau);

    %%  A generation
    A = zeros(m, 1+2*nffs); 
    for j = 1:m
        t = tau(j);
        row = zeros(1, 2*nffs);
        row(1) = 0.5 * (1 - cos(2*pi*t));

        for i = 3:(nffs + 2)
            idx = 2*(i - 3) + 2;
            ca = 0.5*((-1)^i - 1)*cos(pi*t) - ...
                0.5*((-1)^i + 1)*cos(2*pi*t) + cos(i*pi*t);
            cb = i*((-1)^i - 1)*sin(pi*t) - ...
                i*((-1)^i + 1)*sin(2*pi*t) + sin(i*pi*t);
            row(idx)   = ca;
            row(idx+1) = cb;
        end

        A(j,:) = row;
    end

%% Cf generation
    Cf = zeros(m,3);
    for index = 1:3
        % Boundary conditions
        f0 = Xi(index); f1 = Xf(index);
        df0 = Xi(index+3); df1 = Xf(index+3);


        for j = 1:m
            t = tau(j);
            Cf(j,index) = 0.5 * (f0 - f1) * cos(pi*t) + ...
                0.5/pi * (df0 - df1) * sin(pi*t) + ...
                0.5 * (f0 + f1) * cos(pi*t) + ...
                0.25/pi * (df0 + df1) * sin(pi*t);
        end
    end
    % Xx = [a0; reshape(abmatrix', [], 1)];
    % xvec = Ax * Xx + Cf;
end

%% sample run

tou = linspace(0, 1, 100)';
a0 = 1;
coeffs = [0.5, 0.2;
          0.3, 0.3;
          0.2, 0.5]; 
[f, f_dot, f_ddot] = findFFS(tou, a0, coeffs);

[f0, f1, df0, df1] = BCFFS(a0, coeffs);
fprintf("f(0) = %.4f, f(1) = %.4f\n", f0, f1);
fprintf("f'(0) = %.4f, f'(1) = %.4f\n", df0, df1);

[Ax, Cf, xvec, Xx] = buildmatrices(tou, a0, coeffs);
figure;
subplot(3,1,1); plot(tou, f, 'b', tou, xvec, 'r--'); legend('FFS', 'Matrix'); title('f(\tau)');
subplot(3,1,2); plot(tou, f_dot); title('f''(\tau)');
subplot(3,1,3); plot(tou, f_ddot); title('f''''(\tau)');
