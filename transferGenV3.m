function app = transferGenV3
% Creating transferGen programmatic app
fontsize = 16;


% Colors
app.blue = [0.07,0.62,1.00];
app.orange = [0.988, 0.38, 0];
app.purple = [0.659, 0, 1];


% Define all properties.

app.s_data = []; app.t_data = [];
app.s_orbit = []; app.t_orbit = [];
app.CCs = []; app.CCt = [];
app.Ci = []; app.Cf = []; % <--- Display these (J constants)
app.Cx = []; app.Cy = []; %app.Cz = [];
app.dCx = []; app.dCy = []; %app.dCz = [];
app.ddCx = []; app.ddCy = []; %app.ddCz = [];
app.opts113 = [];
app.mu = 1.215058560962404E-2;   % <--- Display this (maybe)?
app.Xstarting = []; app.Xtarget = [];
app.tstarting = []; app.ttarget = [];
app.X0 = []; app.Xf = []; % <---- Save this (no need to display)
app.s_point = []; app.t_point = [];
app.TOF = [];
app.Tp = []; app.Xp = []; app.U = [];
app.m = 100;
app.m = 500;
app.w = [];
app.iPeriod = []; app.fPeriod = []; % <---- periods of initial/final orbits
app.transferArc = [];
app.FouriertransferArc = [];
app.A = []; app.Cx = []; app.Cy = []; app.Cz = [];
app.nffs = 100;
Thrust = [];
deltaV = [];
app.stopRequested = 0;

app.opts113 = odeset('RelTol',1e-10,'AbsTol',1e-11);
app.odefcn = @(~,x) CR3BPMC2D(x,app.mu);

x = [];
y = [];


scalefactor = 1;

% dimensional units

T_em = 2.361 * 10^6;
D = 3.850*10^5;
n = 2*pi/T_em;
M = 1000;

% Initialize FFS algorithm
t = linspace(0,1,app.m)';
[app.A, app.dA, app.ddA] = getAmatrix(t,app.nffs);


fmc_optsfcn = optimoptions("fmincon","MaxFunctionEvaluations", ...
    3000,"Algorithm","interior-point",'OutputFcn',@outputfcn, ...
    'Display','iter-detailed','EnableFeasibilityMode',true);

r1 = zeros(app.m,1);
r2 = zeros(app.m,1);

app.Tmax = 1000;

%% Specifying overall interface structure
app.fig = uifigure;
app.fig.Name = "transferGen";


%app.fig.Position = [0           0        1536         864];
app.fig.WindowStyle = "alwaysontop";
app.fig.Theme = 'light';


app.G1 = uigridlayout(app.fig,[2,2]);
app.G1.ColumnWidth = {'1x','2x'};
% app.G2 = uigridlayout(app.G1,[2,1]);
% app.G3 = uigridlayout

app.P1 = uipanel(app.G1);
app.GP1 = uigridlayout(app.P1,[8,3]);
app.G1.RowHeight = {'2x','1x'};






%% 1st panel: starting/target position selection menu and generate button

app.orbit0InfoLabel = uilabel(app.GP1 ...
    ,'Text','Initial Position Info' ...
    ,'FontSize',fontsize ...
    ,'HorizontalAlignment','center' ...
    ,'FontColor', [0.07,0.62,1.00] ...
    ,'FontWeight','bold');
app.orbit0InfoLabel.Layout.Column = [1,3];

% Input data

%Labels:
% Starting orbit selection label
app.orbit0FamilyLabel = uilabel(app.GP1 ...
    ,'Text', 'Orbit Family:' ...
    ,'FontSize', fontsize ...
    ,'HorizontalAlignment','center');

% Starting orbit
app.orbit0JCLabel = uilabel(app.GP1 ...
    ,'Text', 'Jacobi Constant:' ...
    ,'FontSize', fontsize ...
    ,'HorizontalAlignment','center' ...
    ,'Visible','off');

app.orbit0PosLabel = uilabel(app.GP1 ...
    ,'Text', 'Initial Position Along Orbit:' ...
    ,'FontSize', fontsize ...
    ,'HorizontalAlignment','center' ...
    , 'WordWrap','on' ...
    , 'Visible', 'off');


app.StartingOrbitDropDown = uidropdown(app.GP1,"Placeholder", "<select>" ...
    ,'FontSize',fontsize ...
    ,'ValueChangedFcn',@orbit0DDChanged);
app.StartingOrbitDropDown.Layout.Column  = 1;


app.orbit0JCSlider = uislider(app.GP1,"Visible","off" ...
    ,"ValueChangingFcn",@JacobiConstantSliderValueChanging);


app.orbit0PosSlider = uislider(app.GP1,"MajorTickLabels", ["0","T/4", "T/2", "3T/4","T"] ...
    ,"MajorTicks",[0;0.25;0.5;0.75;1] ...
    ,"Limits",[0,1] ...
    ,'ValueChangingFcn',@pos0changing ...
    ,'Visible','off');


% Target orbit info
% Target orbit info label

app.orbitfInfoLabel = uilabel(app.GP1 ...
    ,'Text','Target Position Info' ...
    ,'FontSize',fontsize ...
    ,'HorizontalAlignment','center' ...
    ,'FontColor', [0.07,0.62,1.00] ...
    ,'FontWeight','bold');
app.orbitfInfoLabel.Layout.Column = [1,3];
app.orbitfInfoLabel.Layout.Row = 4;


% Target Orbit Family Selection Label
app.orbitfFamilyLabel = uilabel(app.GP1 ...
    ,'Text', 'Orbit Family:' ...
    ,'FontSize', fontsize ...
    ,'HorizontalAlignment','center');
app.orbitfFamilyLabel.Layout.Row = 5;

% Target Orbit Jacobi Constant Selection Label
app.orbitfJCLabel = uilabel(app.GP1 ...
    ,'Text', 'Jacobi Constant:' ...
    ,'FontSize', fontsize ...
    ,'HorizontalAlignment','center' ...
    ,'Visible','off');

% Target Position Location Label
app.orbitfPosLabel = uilabel(app.GP1 ...
    ,'Text', 'Initial Position Along Orbit:' ...
    ,'FontSize', fontsize ...
    ,'HorizontalAlignment','center' ...
    , 'WordWrap','on' ...
    ,'Visible','off');



% Target orbit selection drop down object
app.TargetOrbitDropDown = uidropdown(app.GP1,"Placeholder", "<select>" ...
    ,'FontSize',fontsize ...
    ,'ValueChangedFcn',@orbitfDDChanged);
app.TargetOrbitDropDown.Layout.Column  = 1;


% Target Orbit Jacobi Constant Selection Slider Object
app.orbitfJCSlider = uislider(app.GP1 ...
    ,"ValueChangingFcn",@TargetJacobiConstantSliderValueChanging ...
    ,'Visible','off');


% Target position location slider object
app.orbitfPosSlider = uislider(app.GP1,"MajorTickLabels", ["0","T/4", "T/2", "3T/4","T"] ...
    ,"MajorTicks",[0;0.25;0.5;0.75;1] ...
    ,"Limits",[0,1] ...
    ,'ValueChangingFcn',@posfchanging ...
    ,'Visible','off');

% Generate transfer button object
app.generate = uibutton(app.GP1 ...
    ,"Text",'Generate Transfer' ...
    ,'ButtonPushedFcn',@generateButtonPushed);
app.generate.Layout.Column = 2;
app.generate.Layout.Row = 8;


% Generate transfer label
app.generateLabel = uilabel(app.GP1 ...
    ,'Text','Transfer Generation' ...
    ,'FontSize',fontsize ...
    ,'HorizontalAlignment','center' ...
    ,'FontColor', [0.988, 0.38, 0] ...
    ,'FontWeight','bold');
app.generateLabel.Layout.Column = [1,3];
app.generateLabel.Layout.Row = 7;


app.stopbutton = uibutton(app.GP1 ...
    ,'Text', 'Cancel' ...
    ,'FontSize',fontsize ...
    ,'ButtonPushedFcn', @stopExecution);






% Initializing Axis For Plotting
app.UIAxes = uiaxes(app.G1);
app.UIAxes.Layout.Row = [1 2];


app.UIAxes2 = uiaxes('Parent',app.G1);
app.UIAxes2.Visible = 'off';


% Window  Settings
app.fig.WindowState = "maximized";
app.fig.WindowStyle = "alwaysontop";

%% Panel 2 --- Will display critical information:
% - Jacobi Constants of starting and target orbits
% - Indices of starting and target orbits in their respective .csv
% files (and potentially the ability to precisely enter the required index)
% - Required delta V for maneuver
% - Maneuver time of flight
% - Nondimensional departure and arrival times relative to ICs on .csv
% files.



app.P2 = uipanel(app.G1);

app.GP2 = uigridlayout(app.P2,[8,4]);
app.G2.RowHeight = {'2x','1x'};

app.P2.Layout.Row = 2;
app.P2.Layout.Column = 1;

%Orbit Headers
app.generateLabel = uilabel(app.GP2 ...
    ,'Text','Orbit 1' ...
    ,'FontSize',fontsize ...
    ,'HorizontalAlignment','center' ...
    ,'FontColor', [0,0, 0] ...
    ,'FontWeight','bold');
app.generateLabel.Layout.Column = [1,2];
app.generateLabel.Layout.Row = 1;

app.generateLabel = uilabel(app.GP2 ...
    ,'Text','Orbit 2' ...
    ,'FontSize',fontsize ...
    ,'HorizontalAlignment','center' ...
    ,'FontColor', [0,0, 0] ...
    ,'FontWeight','bold');
app.generateLabel.Layout.Column = [3,4];
app.generateLabel.Layout.Row = 1;


% Jacobi Constants and Indices

app.CiLabel = uilabel(app.GP2, 'Text', 'Jacobi Constant:');
app.CiLabel.Layout.Column = [1,2];
app.CiLabel.Layout.Row = 2;

app.IndexiLabel = uilabel(app.GP2, 'Text', 'Index of Orbit:');
app.IndexiLabel.Layout.Column = [1,2];
app.IndexiLabel.Layout.Row = 3;

app.CfLabel = uilabel(app.GP2, 'Text', 'Jacobi Constant:');
app.CfLabel.Layout.Column = [3,4];
app.CfLabel.Layout.Row = 2;

app.IndexfLabel = uilabel(app.GP2, 'Text', 'Index of Orbit:');
app.IndexfLabel.Layout.Column = [3,4];
app.IndexfLabel.Layout.Row = 3;

% Required dV

app.generateLabel = uilabel(app.GP2 ...
    ,'Text','Tranfer Characteristics' ...
    ,'FontSize',fontsize ...
    ,'HorizontalAlignment','center' ...
    ,'FontColor', [0,0, 0] ...
    ,'FontWeight','bold');
app.generateLabel.Layout.Column = [2,3];
app.generateLabel.Layout.Row = 4;

app.dVLabel = uilabel(app.GP2, 'Text', 'Transfer dV:');
app.dVLabel.Layout.Column = [2,3];
app.dVLabel.Layout.Row = 5;

app.TOFLabel = uilabel(app.GP2, 'Text', 'Transfer Time:');
app.TOFLabel.Layout.Column = [2,3];
app.TOFLabel.Layout.Row = 6;



%% Extracting list of files from specified folder.

if ismac
    addpath("OrbitDataEarthMoon/")
    files = dir("OrbitDataEarthMoon/*.csv");
    list_of_families = ({files.name}');
    %list_of_families = cellstr(list_of_families)';
    array_of_families = cell(1,size(list_of_families,1));
    for i = 1:size(list_of_families,1)
        array_of_families{i} = char(list_of_families(i,:));
    end
elseif isunix
    addpath("OrbitDataEarthMoon\")
    list_of_families = ls("OrbitDataEarthMoon\*.csv");
elseif ispc
    addpath("OrbitDataEarthMoon\")
    list_of_families = ls("OrbitDataEarthMoon\*.csv");
    array_of_families = cell(1,size(list_of_families,1));
    for i = 1:size(list_of_families,1)
        array_of_families{i} = list_of_families(i,:);
    end
else
    disp('Platform not supported')
end

app.StartingOrbitDropDown.Items = array_of_families;
app.TargetOrbitDropDown.Items = array_of_families;


%% plotting Earth, Moon and Lagrange Points
Re = 6378/384400;
Rm = 1737/389703;

% [Xe,Ye,Ze] = ellipsoid(app.UIAxes,-app.mu,0,0,Re,Re,Re,400);
% [Xm,Ym,Zm] = ellipsoid(app.UIAxes,1-app.mu,0,0,Rm,Rm,Rm,400);

xL1 = Lptpos(app.mu,1);
xL2 = Lptpos(app.mu,2);
xL3 = Lptpos(app.mu,3);

hold (app.UIAxes, 'on')
draw_shaded_circle(app.UIAxes,[-app.mu,0,0],Re, app.blue)
draw_shaded_circle(app.UIAxes,[1-app.mu,0,0],Rm, [0.1,0.1,0.1])
plot(app.UIAxes, xL1, 0, 'rx')
plot(app.UIAxes, xL2, 0, 'rx')
plot(app.UIAxes, xL3, 0, 'rx')
plot(app.UIAxes,1/2,sqrt(3)/2,'rx','LineWidth',1)
plot(app.UIAxes,1/2,-sqrt(3)/2,'rx','LineWidth',1)
hold (app.UIAxes,'off')


xlabel(app.UIAxes,'$x$ [DU]','FontSize',14,'Interpreter','latex')
ylabel(app.UIAxes,'$y$ [DU]','FontSize',14,'Interpreter','latex')
% zlabel(app.UIAxes,'$z$ [DU]','FontSize',14,'Interpreter','latex')
title(app.UIAxes,'')
axis(app.UIAxes,'equal')
%view(app.UIAxes,[0,0,1])
grid(app.UIAxes,"on")
box(app.UIAxes,"on")





%% Callback functions

% Selecting new starting orbit family
    function orbit0DDChanged(src, event)
        addpath("OrbitDataEarthMoon\")

        value = event.Value;
        if ~app.orbit0JCLabel.Visible
            app.orbit0JCLabel.Visible = "on";
            app.orbit0JCSlider.Visible = 'on';
        end

        s_data_name = value;
        app.s_data = readmatrix(s_data_name);
        app.CCs = app.s_data(:,8);
        Cmin = min(app.CCs);
        Cmax = max(app.CCs);
        app.orbit0JCSlider.Limits=[Cmin,Cmax];
        %app.EditField.Value = Cmin;

    end

% Selecting new target orbit family
    function orbitfDDChanged(src, event)
        addpath("OrbitDataEarthMoon\")

        value = event.Value;

        if ~app.orbitfJCLabel.Visible
            app.orbitfJCLabel.Visible = 'on';
            app.orbitfJCSlider.Visible = 'on';
        end

        t_data_name = value;
        app.t_data = readmatrix(t_data_name);
        app.CCt = app.t_data(:,8);
        Cmin = min(app.CCt);
        Cmax = max(app.CCt);
        app.orbitfJCSlider.Limits=[Cmin,Cmax];

        %app.EditField.Value = Cmin;


    end

% Selecting new starting Jacobi Constant
    function JacobiConstantSliderValueChanging(src, event)

        changingValue = event.Value;
        value = changingValue;

        if ~app.orbit0PosLabel.Visible
            app.orbit0PosLabel.Visible = "on";
            app.orbit0PosSlider.Visible = "on";
        end

        [~,ind] = min(abs(value - app.CCs));
        C = app.CCs(ind);
        indexincsv = app.s_data(ind,1);
        app.CiLabel.Text = sprintf('Jacobi Constant: %.9f', C);
        app.IndexiLabel.Text = sprintf('Index of Orbit: %d', indexincsv);
        %app.EditField.Value = C;
        X0 = app.s_data(ind,[2,3,5,6])';
        T = app.s_data(ind,9);
        app.iPeriod = T;
        %tic;
        [app.tstarting,app.Xstarting] = ode113(app.odefcn,[0,T],X0,app.opts113);
        %toc;
        x = app.Xstarting(:,1); y = app.Xstarting(:,2);
        delete(app.s_orbit);
        delete(app.s_point);
        delete(app.transferArc)
        hold(app.UIAxes,'on')
        app.s_orbit = plot(app.UIAxes,x,y,'Color',[0.07,0.62,1.00],'LineWidth',2);
        hold(app.UIAxes,'off')

    end

% Selecting new target Jacobi Constant
    function TargetJacobiConstantSliderValueChanging(src, event)

        changingValue = event.Value;
        value = changingValue;

        if ~app.orbitfPosLabel.Visible
            app.orbitfPosLabel.Visible = "on";
            app.orbitfPosSlider.Visible = "on";
        end

        [~,ind] = min(abs(value - app.CCt));
        C = app.CCt(ind);
        indexincsv = app.t_data(ind,1);
        % app.EditField.Value = C;
        app.CfLabel.Text = sprintf('Jacobi Constant: %.9f', C);
        app.IndexfLabel.Text = sprintf('Index of Orbit: %d', indexincsv);

        X0 = app.t_data(ind,[2,3,5,6])';
        T = app.t_data(ind,9);
        app.fPeriod = T;
        %tic;
        [app.ttarget,app.Xtarget] = ode113(app.odefcn,[0,T],X0,app.opts113);
        %toc;
        x = app.Xtarget(:,1); y = app.Xtarget(:,2);
        delete(app.t_orbit);
        delete(app.t_point);
        delete(app.transferArc)
        hold(app.UIAxes,'on')
        app.t_orbit = plot(app.UIAxes,x,y,'Color',[0.659, 0, 1],'LineWidth',2);
        hold(app.UIAxes,'off')
    end

% Selecting new starting position location
    function pos0changing(src,event)
        changingValue = event.Value;

        if ~isempty(app.tstarting)
            tau = app.tstarting/app.tstarting(end);
            [~,ind] = min(abs(changingValue-tau));
            app.X0 = app.Xstarting(ind,:);

            delete(app.s_point);
            delete(app.transferArc)
            hold(app.UIAxes,"on")
            app.s_point = plot(app.UIAxes, app.X0(1),app.X0(2),'go','LineWidth',3,'MarkerFaceColor','g');
            hold(app.UIAxes,"off")

        end

    end

% Selecting new target position location
    function posfchanging(src,event)
        changingValue = event.Value;

        if ~isempty(app.ttarget)
            tau = app.ttarget/app.ttarget(end);
            [~,ind] = min(abs(changingValue-tau));
            app.Xf = app.Xtarget(ind,:);

            delete(app.t_point);
            delete(app.transferArc)
            hold(app.UIAxes,"on")
            app.t_point = plot(app.UIAxes, app.Xf(1),app.Xf(2),'ro','LineWidth',3,'MarkerFaceColor','r');
            hold(app.UIAxes,"off")

        end

    end

% Generate optimal transfer
    function generateButtonPushed(src,event)

        if ~isempty(app.s_point) & ~isempty(app.t_point)
            % Discretization points
            delete(app.transferArc)
            delete(app.FouriertransferArc)
        
            app.UIAxes.Layout.Row = 1;
            drawnow;
            
            app.UIAxes2.Visible = 'on';
            app.UIAxes2.Layout.Row = 2;
            app.UIAxes2.Layout.Column = 2;
            drawnow;

            % Manifold generation ---------------------------------------

            % Generate initial conditions for var eq
            I = eye(4);
            Qi = [reshape(I,16,1); app.X0'];
            Qf = [reshape(I,16,1); app.Xf'];

            varfun = @(t,u) var2D(t,u,app.mu); 



            [ti, PHIi] = ode113(varfun,[0,app.iPeriod],Qi,app.opts113);
            [tf, PHIf] = ode113(varfun,[0,app.fPeriod],Qf,app.opts113);

            XXi = PHIi(:,17:20);
            XXf = PHIf(:,17:20);

            % phi_i = reshape(PHIi(end,1:16),[4,4]);
            % phi_f = reshape(PHIf(end,1:16),[4,4]);
            % 
            % [Vi,Lambdai] = eigs(phi_i);
            % [Vf, Lambdaf] = eigs(phi_f);
            % 
            % 
            % [~,index] = max(real(diag(Lambdai)));
            % 
            % u_unst = Vi(:,index);
            % 
            % [~,index] = min(real(diag(Lambdaf)));
            % 
            % u_stab = Vf(:,index);
            % 
            % W_unst_0 = app.X0' + 5e-6*u_unst;
            % W_stab_0 = app.Xf' - 5e-6*u_stab;
            % 
            % [~,W_unst] = ode113(app.odefcn, [0,7*app.iPeriod],W_unst_0,app.opts113);
            % [~,W_stab] = ode113(app.odefcn, [0,-7*app.fPeriod],W_stab_0,app.opts113);
            % 
            % hold(app.UIAxes,'on')
            % 
            % plot(app.UIAxes, W_unst(:,1), W_unst(:,2),'g','LineWidth',2)
            % plot(app.UIAxes, W_stab(:,1), W_stab(:,2),'r','LineWidth',2)    

            % Manifold generation ---------------------------------------


            orbit1 = interp1(ti/app.iPeriod,XXi,t);
            orbit2 = interp1(tf/app.fPeriod, XXf,t);

            %
            [init_guess,app.w] = interp_arc(orbit1,orbit2,app.m);
            x = init_guess(:,1);
            y = init_guess(:,2);


            % We will start with this
            findBestTf(1);

            app.qx = pinv(app.A) * (x - app.Cx);
            app.qy = pinv(app.A) * (y - app.Cy);

            % Time is NOT a design variable
            % z0 = [app.qx;app.qy];

            % Time is a design variable:
            z0 = [app.qx; app.qy; app.TOF * scalefactor];

            % Linear constraint on maneuver time:
            e = eye(length(z0));
            e = e(end,:);

            % no TOF
            % problem = struct('objective',@wrappedCostFcn,'x0',z0 ...
            %     ,'Aineq',[],'bineq', [], 'Aeq', [], 'beq', [] ...
            %     ,'lb', [], 'ub', [], 'nonlcon', @(z) confunc(z,0.01866,0.00478) ...
            %     ,'solver','fmincon', 'options', fmc_optsfcn);

            % with nonlinear constraint
            problem = struct('objective',@wrappedCostFcn,'x0',z0 ...
                ,'Aineq',e,'bineq', app.Tmax, 'Aeq', [], 'beq', [] ...
                ,'lb', [], 'ub', [], 'nonlcon', @(z) confunc(z,0.01866,0.00478) ...
            ,'solver','fmincon', 'options', fmc_optsfcn);

            app.stopRequested = 0;
            [z_opt] = fmincon(problem);

            z = z_opt;
            app.TOF = z(end);
            N = 2*app.nffs+1;
            app.qx = z(1:N);
            app.qy = z(N+1:end-1);
  

            x = app.A * app.qx + app.Cx;
            y = app.A * app.qy + app.Cy;
    


            hold(app.UIAxes,'on')
            app.transferArc = plot(app.UIAxes,x,y,'LineWidth',2, ...
                'Color',app.orange);
            % app.FouriertransferArc = plot3(app.UIAxes, xtest,ytest,ztest ...
            %     ,'r--','LineWidth',2);
            hold(app.UIAxes,'off')
        end

    end

% Cost function Definition
    function J = transferGenCostFcn(Z)
        % z:  3(2nffs-3) x 1
        % = [TOF   qx'   qy'  qz']'
        % where qx = [a_x0 a_x3 b_x3 a_x4 b_x4 ... a_xN, b_xN]'
        %       qy = [a_y0 a_y3 b_y3 a_y4 b_y4 ... a_yN, b_yN]'
        %       qz = [a_z0 a_z3 b_z3 a_z4 b_z4 ... a_zN, b_zN]'

        app.TOF = Z(end) / scalefactor;
        %app.TOFLabel.Text = sprintf('Transfer Time: %f', app.TOF);

        N = 2*app.nffs+1;
        app.qx = Z(1:N);
        app.qy = Z(N+1:2*N);

       [Cf,dCf, ddCf] = findCf2D(t,app.X0, app.Xf, app.m);
       app.Cx = Cf(:,1); app.Cy = Cf(:,2);
       app.dCx = dCf(:,1); app.dCy = dCf(:,2); 
       app.ddCx = ddCf(:,1); app.ddCy = ddCf(:,2);


        x = app.A * app.qx + app.Cx;
        y = app.A * app.qy + app.Cy;


        xd = 1/app.TOF * (app.dA * app.qx + app.dCx);
        yd = 1/app.TOF * (app.dA * app.qy + app.dCy);
        %zd = 1/app.TOF * (app.dA * app.qz + app.dCz);


        xdd = 1/app.TOF^2 * (app.ddA * app.qx + app.ddCx);
        ydd = 1/app.TOF^2 * (app.ddA * app.qy + app.ddCy);



        r1 = sqrt((x+app.mu).^2 +y.^2);
        r2 = sqrt((x-(1-app.mu)).^2 +y.^2);

        ux = xdd - 2*yd - x +...
            (1-app.mu)*(x+app.mu)./r1.^3 + ...
            app.mu*(x-1+app.mu)./r2.^3;

        uy = ydd + 2*xd - y +...
            (1-app.mu)*y./r1.^3 + app.mu*y./r2.^3  ;



        for k = 1:app.m
            U = [ux(k);uy(k)];
            L(k) = U'*U;
            Thrust(k) = sqrt(L(k))*M*n^2*D;
        end

        J = trapz(t*app.TOF, L);


        deltaV =  1/M *trapz(t*app.TOF*T_em/(2*pi),Thrust);

        sqrt(J);

        %sqrt(J)

        J = deltaV;
    end


    function [c,ceq] = confunc(Z, re, rm)

        app.TOF = Z(end) / scalefactor;
        

        N = 2*app.nffs+1;
        app.qx = Z(1:N);
        app.qy = Z(N+1:2*N);

       [Cf,dCf, ddCf] = findCf2D(t,app.X0, app.Xf, app.m);
       app.Cx = Cf(:,1); app.Cy = Cf(:,2);
       app.dCx = dCf(:,1); app.dCy = dCf(:,2); 
       app.ddCx = ddCf(:,1); app.ddCy = ddCf(:,2);


        x = app.A * app.qx + app.Cx;
        y = app.A * app.qy + app.Cy;


        xd = 1/app.TOF * (app.dA * app.qx + app.dCx);
        yd = 1/app.TOF * (app.dA * app.qy + app.dCy);
        %zd = 1/app.TOF * (app.dA * app.qz + app.dCz);


        xdd = 1/app.TOF^2 * (app.ddA * app.qx + app.ddCx);
        ydd = 1/app.TOF^2 * (app.ddA * app.qy + app.ddCy);



        r1 = sqrt((x+app.mu).^2 +y.^2);
        r2 = sqrt((x-(1-app.mu)).^2 +y.^2);

        ux = xdd - 2*yd - x +...
            (1-app.mu)*(x+app.mu)./r1.^3 + ...
            app.mu*(x-1+app.mu)./r2.^3;

        uy = ydd + 2*xd - y +...
            (1-app.mu)*y./r1.^3 + app.mu*y./r2.^3  ;



        for k = 1:app.m
            U = [ux(k);uy(k)];
            L(k) = U'*U;
            Thrust(k) = sqrt(L(k))*M*n^2*D;
        end
        c1 = re - r1;
        c2 = rm - r2;
        c3 = Thrust'*1000 - 1;
        c = [];
        ceq = [];
    end

    function stop = outputfcn(Z,optimValues,state)
        stop = false;

        if any(strcmp(state, {'iter','interrupt','done','init'}))

            deltaV_label = ['Transfer dV: ', num2str(deltaV*1000) ' m/s'];
            app.dVLabel.Text = deltaV_label;
            TOF_days = app.TOF*T_em/(2*pi*3600*24);

            time_label = ['Transfer Time:', num2str(TOF_days) ' days'];
            app.TOFLabel.Text = time_label;

            delete(app.transferArc)


            hold(app.UIAxes, 'on');
            app.transferArc = plot(app.UIAxes, x, y, ...
                'LineWidth', 2, 'Color', app.orange);

            hold(app.UIAxes, 'off');


            % Plotting Thrust vs Time Plot in lower right corner
            %app.UIAxes2.Layout.Row = 2;
            app.UIAxes2.Layout.Column = 2;

            plot(app.UIAxes2,t, Thrust*1000,'LineWidth',2)
            grid(app.UIAxes2,'on')
            xlabel(app.UIAxes2,'$\tau$','FontSize',fontsize ...
                ,'Interpreter','latex')
            ylabel(app.UIAxes2,'Thrust $\mathcal{T}$ [N]' ...
                ,'FontSize',fontsize,'Interpreter','latex')
            drawnow expose;


        end
    end

    function stopExecution(src,event)
        app.stopRequested = true;
    end

    function J = wrappedCostFcn(Z)
        J = transferGenCostFcn(Z);

        deltaV_label = ['Transfer dV: ', num2str(deltaV*1000) ' m/s'];
        app.dVLabel.Text = deltaV_label;

        app.TOF = Z(end) / scalefactor;
        TOF_days = app.TOF*T_em/(2*pi*3600*24);

        time_label = ['Transfer Time:', num2str(TOF_days) ' days'];
        app.TOFLabel.Text = time_label;

        delete(app.transferArc)


        hold(app.UIAxes, 'on');
        app.transferArc = plot(app.UIAxes, x, y, ...
            'LineWidth', 2, 'Color', app.orange);

        hold(app.UIAxes, 'off');


        % Plotting Thrust vs Time Plot in lower right corner
        %app.UIAxes2.Layout.Row = 2;
        app.UIAxes2.Layout.Column = 2;

        plot(app.UIAxes2,t, Thrust*1000,'LineWidth',2)
        grid(app.UIAxes2,'on')
        xlabel(app.UIAxes2,'$\tau$','FontSize',fontsize ...
            ,'Interpreter','latex')
        ylabel(app.UIAxes2,'Thrust $\mathcal{T}$ [N]' ...
            ,'FontSize',fontsize,'Interpreter','latex')
        drawnow expose;


        if app.stopRequested
            error('Optimization interrupted by user.');
        end

    end

    function findBestTf(pl)
        % plot - 0 or 1, shows plot of J0 vs Tf.
        %Z0 = [app.qx; app.qy; app.qz];
        Tfrange = linspace(0.1,app.Tmax,50);

        for k = 1:length(Tfrange)
            app.TOF = Tfrange(k);

            [Cf,dCf, ddCf] = findCf2D(t,app.X0, app.Xf, app.m);
            app.Cx = Cf(:,1); app.Cy = Cf(:,2);
            app.dCx = dCf(:,1); app.dCy = dCf(:,2); 
            app.ddCx = ddCf(:,1); app.ddCy = ddCf(:,2);

            app.qx = pinv(app.A) * (x - app.Cx);
            app.qy = pinv(app.A) * (y - app.Cy);

            Z0 = [app.qx; app.qy; app.TOF * scalefactor];
            [~] = transferGenCostFcn(Z0);
            dV(k) = deltaV;
        end

        [min_dV,ind] = min(dV);
        app.TOF = Tfrange(ind);
        deltaV = min_dV;

        [Cf,dCf, ddCf] = findCf2D(t,app.X0, app.Xf, app.m);
        app.Cx = Cf(:,1); app.Cy = Cf(:,2);
        app.dCx = dCf(:,1); app.dCy = dCf(:,2); 
        app.ddCx = ddCf(:,1); app.ddCy = ddCf(:,2);


        deltaV_label = ['Transfer dV: ', num2str(deltaV*1000) ' m/s'];
        app.dVLabel.Text = deltaV_label;
        TOF_days = app.TOF*T_em/(2*pi*3600*24);

        time_label = ['Transfer Time:', num2str(TOF_days) ' days'];
        app.TOFLabel.Text = time_label;

        delete(app.transferArc)
        hold(app.UIAxes,'on')
        app.transferArc = plot(app.UIAxes,x,y,'LineWidth',2, ...
            'Color',app.orange);
        hold(app.UIAxes,'off')

        if pl
            figure();
            hold on
            plot(Tfrange,dV,'LineWidth',2)
            plot(Tfrange(ind),dV(ind),'ro')
            hold off
            xlabel ('Tf')
            ylabel ('J')
            grid on
            drawnow;
        end

    end

    function draw_shaded_circle(ax, center, radius, color)
        % draw_shaded_circle(center, radius, color)
        % Draws a filled circle at 'center' = [x, y] with given 'radius'
        % and optional 'color' (default is light blue).
        
        if nargin < 4
            color = [0.3, 0.6, 1]; % default color: light blue
        end
    
        theta = linspace(0, 2*pi, 100);  % 100 points around the circle
        u = center(1) + radius * cos(theta);
        v = center(2) + radius * sin(theta);
    
        fill(ax, u, v, color, 'EdgeColor', 'none');  % filled, no border
        %axis equal;  % keep aspect ratio correct
    end

function [Cf,Cfd,Cfdd] = findCf2D(tou, Xi, Xf, m)
    
    Cf = zeros(m,2);
    Cfd = zeros(m,2);
    Cfdd = zeros(m,2);
    for index = 1:2
        % Boundary conditions
        f0 = Xi(index); f1 = Xf(index);
        df0 = app.TOF * Xi(index+2); df1 =  app.TOF * Xf(index+2);


        for j = 1:m
            tt = tou(j);
            Cf(j,index) = 0.5 * (f0 - f1) * cos(pi*tt) + ...
                0.5/pi * (df0 - df1) * sin(pi*tt) + ...
                0.5 * (f0 + f1) * cos(2*pi*tt) + ...
                0.25/pi * (df0 + df1) * sin(2*pi*tt);
            Cfd(j,index) = -(pi/2) * (f0-f1) * sin(pi*tt) + ...
            (df0-df1) * cos(pi*tt) * 0.5 -pi*(f0+f1) * sin(2*pi*tt) + ...
            (df0+df1)*cos(2*pi*tt) * 0.5;

            Cfdd(j,index) = -(pi*pi*0.5)*(f0-f1)*cos(pi*tt) - ...
            0.5 * pi * (df0-df1) * sin(pi*tt) - ...
            2*  pi * pi * (f0+f1) * cos(2*pi*tt) - ...
            pi * (df0+df1) * sin(2*pi*tt);

        end
    end
end

end


