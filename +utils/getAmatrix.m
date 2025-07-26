function [AAA] = getAmatrix(m, n)
    % after merge
    tau = linspace(0,1,m);

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
    AAA = [A;dA;ddA];

end
