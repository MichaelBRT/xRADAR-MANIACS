clear, clc, close all





folder_name = 'PlanarOrbitData';
folder_path = fullfile(pwd,folder_name);


files = dir(folder_path);
files = files(~[files.isdir]);  % Remove '.' and '..'

N = numel(files);

for i = 1:N
    orbit_name = files(i).name;
    output_file = ['filtered PlanarOrbitData\', orbit_name];
    if ~exist(output_file,"file")

        % Pull orbit data
        orbit_file = ['PlanarOrbitData\', orbit_name];
        data = readmatrix(orbit_file);

        mu = data(1,end);

        % We want only Jacobi Constants C < C_L3
        xL3 = Lptpos(mu,3);
        C_L3 = Jconst([xL3;0;0;0]);

        idx = find(data(:,8) > C_L3);
        data = data(idx,:);

        % Now we wish to replace the JPL stability index with 1/2*(lambda + 1/lambda)
        if ~isempty(data)
            vareqn = @(t,x) var2D(t,x,mu);
            % Smallest possible tolerance
            varopt = odeset('RelTol',3e-14,'AbsTol',1e-14);
            phi_0 = reshape(eye(4),[16,1]);

            for i = 1:size(data,1)
                X0 = data(i,[2,3,5,6])';
                T = data(i,9);
                Y0 = [phi_0;X0];
                [~,Yspan] = ode113(vareqn,[0,T],Y0,varopt);
                M = reshape(Yspan(end,1:16)',[4,4]);
                data(i,end-1) = stability_index(M);
            end

            % Now we wish to extract only the unstable orbits
            idx = find(abs(data(:,end-1))>1);
            data = data(idx,:);

            % Now we wish to save the orbit to a new file under filtered
            % PlanarOrbitData
            if ~isempty(data)
                writematrix(data, output_file);
            end
        end
    end

end




function C = Jconst(X)
% Function calculates Jacobi constant given a state
x = X(1);
y = X(2);
xdot = X(3);
ydot = X(4);

mu = 0.012150584270572;
mu1=1-mu;
mu2=  mu;

r1 = ((x+mu2)^2 + y^2)^(1/2);
r2 = ((x-mu1)^2 + y^2)^(1/2);

U = -1/2*(x^2+y^2)-mu1/r1-mu2/r2; %-1/2*mu1*mu2;

C = -(xdot^2 + ydot^2)-2*U;
end

function I = stability_index(M)
V  = eigs(M);
[~,inds] = maxk(real(V)-1,2,'ComparisonMethod','abs');
I = 1/2*sum(V(inds));
end