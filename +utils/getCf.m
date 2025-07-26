function [CCC] = getCf(tou, Xi, Xf, TOM)
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
    CCC = [Cf,Cfd,Cfdd];
 end