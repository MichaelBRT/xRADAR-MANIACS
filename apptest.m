function app = transferGenV2
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
    app.Cx = []; app.Cy = []; app.Cz = [];
    app.dCx = []; app.dCy = []; app.dCz = [];
    app.ddCx = []; app.ddCy = []; app.ddCz = [];
    app.opts113 = [];
    app.mu = 1.215058560962404E-2;   % <--- Display this (maybe)?
    app.Xstarting = []; app.Xtarget = [];
    app.tstarting = []; app.ttarget = [];
    app.X0 = []; app.Xf = []; % <---- Save this (no need to display)
    app.s_point = []; app.t_point = []; 
    app.TOF = [];
    app.Tp = []; app.Xp = []; app.U = [];
    app.m = 1000;
    app.w = [];
    app.iPeriod = []; app.fPeriod = []; % <---- periods of initial/final orbits
    app.transferArc = [];
    app.FouriertransferArc = [];
    app.A = []; app.Cx = []; app.Cy = []; app.Cz = [];
    app.nffs = 100;
   
    app.opts113 = odeset('RelTol',1e-10,'AbsTol',1e-11);
    app.odefcn = @(~,x) CR3BPMC(x,app.mu);


    % Initialize FFS algorithm
    t = linspace(0,1,app.m)';
    [app.A, app.dA, app.ddA] = getAmatrix(t,app.nffs);


    fmc_optsfcn = optimoptions("fmincon","MaxFunctionEvaluations",100,"Algorithm","sqp");

    r1 = zeros(app.m,1);
    r2 = zeros(app.m,1);

    app.Tmax = 100;
    
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




   

    % Initializing Axis For Plotting
    app.UIAxes = uiaxes(app.G1);
    app.UIAxes.Layout.Row = [1 2];

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
    Re = 6378/389703;
    Rm = 1737/389703;

    [Xe,Ye,Ze] = ellipsoid(app.UIAxes,-app.mu,0,0,Re,Re,Re,400);
    [Xm,Ym,Zm] = ellipsoid(app.UIAxes,1-app.mu,0,0,Rm,Rm,Rm,400);
    
    xL1 = Lptpos(app.mu,1);
    xL2 = Lptpos(app.mu,2);
    xL3 = Lptpos(app.mu,3);

    hold (app.UIAxes, 'on')
    surf(app.UIAxes,Xe,Ye,Ze,'EdgeColor','b','FaceColor','flat')
    app.moon_plot = surf(app.UIAxes,Xm,Ym,Zm,'EdgeColor',[0.1,0.1,0.1],'FaceColor','flat');
    plot3(app.UIAxes, xL1, 0, 0, 'rx')
    plot3(app.UIAxes, xL2, 0, 0, 'rx')
    plot3(app.UIAxes, xL3, 0, 0, 'rx')
    plot3(app.UIAxes,1/2,sqrt(3)/2,0,'rx','LineWidth',1)
    plot3(app.UIAxes,1/2,-sqrt(3)/2,0,'rx','LineWidth',1)
    hold (app.UIAxes,'off')


    xlabel(app.UIAxes,'$x$ [DU]','FontSize',14,'Interpreter','latex')
    ylabel(app.UIAxes,'$y$ [DU]','FontSize',14,'Interpreter','latex')
    zlabel(app.UIAxes,'$z$ [DU]','FontSize',14,'Interpreter','latex')
    title(app.UIAxes,'')
    axis(app.UIAxes,'equal')
    view(app.UIAxes,[0,0,1])
    grid(app.UIAxes,"on")
    box(app.UIAxes,"on")
    pbaspect(app.UIAxes, [1,1,1]);

    

    
    %% Callback functions

% Selecting new starting orbit family
function orbit0DDChanged(src, event)


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
            app.CiLabel.Text = sprintf('Jacobi Constant: %.9f', C);
            app.IndexiLabel.Text = sprintf('Index of Orbit: %d', ind);
            %app.EditField.Value = C;

            X0 = app.s_data(ind,2:7)';
            T = app.s_data(ind,9);
            app.iPeriod = T;
            %tic;
            [app.tstarting,app.Xstarting] = ode113(app.odefcn,[0,T],X0,app.opts113);
            %toc;
            x = app.Xstarting(:,1); y = app.Xstarting(:,2); z = app.Xstarting(:,3);
            delete(app.s_orbit);
            delete(app.s_point);
            delete(app.transferArc)
            hold(app.UIAxes,'on')
            app.s_orbit = plot3(app.UIAxes,x,y,z,'Color',[0.07,0.62,1.00],'LineWidth',2);
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
            % app.EditField.Value = C;
            app.CfLabel.Text = sprintf('Jacobi Constant: %.9f', C);
            app.IndexfLabel.Text = sprintf('Index of Orbit: %d', ind);

            X0 = app.t_data(ind,2:7)';
            T = app.t_data(ind,9);
            app.fPeriod = T;
            %tic;
            [app.ttarget,app.Xtarget] = ode113(app.odefcn,[0,T],X0,app.opts113);
            %toc;
            x = app.Xtarget(:,1); y = app.Xtarget(:,2); z = app.Xtarget(:,3);
            delete(app.t_orbit);
            delete(app.t_point);
            delete(app.transferArc)
            hold(app.UIAxes,'on')
            app.t_orbit = plot3(app.UIAxes,x,y,z,'Color',[0.659, 0, 1],'LineWidth',2);
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
                app.s_point = plot3(app.UIAxes, app.X0(1),app.X0(2),app.X0(3),'go','LineWidth',3,'MarkerFaceColor','g');
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
                app.t_point = plot3(app.UIAxes, app.Xf(1),app.Xf(2),app.Xf(3),'ro','LineWidth',3,'MarkerFaceColor','r');
                hold(app.UIAxes,"off")

            end

    end

    % Generate optimal transfer
    function generateButtonPushed(src,event)
        
        if ~isempty(app.s_point) & ~isempty(app.t_point)
            % Discretization points
            delete(app.transferArc)
            delete(app.FouriertransferArc)

            

            [ti, XXi] = ode113(app.odefcn,[0,app.iPeriod],app.X0,app.opts113);
            [tf, XXf] = ode113(app.odefcn,[0,app.fPeriod],app.Xf,app.opts113);
            orbit1 = interp1(ti/app.iPeriod,XXi,t);
            orbit2 = interp1(tf/app.fPeriod, XXf,t);
    
            %
            [init_guess,app.w] = interp_arc(orbit1,orbit2,app.m);
            x = init_guess(:,1);
            y = init_guess(:,2);
            z = init_guess(:,3);


            %

            [Cf,dCf, ddCf] = findCf(t,app.X0, app.Xf, app.m);
            app.Cx = Cf(:,1); app.Cy = Cf(:,2); app.Cz = Cf(:,3);
            app.dCx = dCf(:,1); app.dCy = dCf(:,2); app.dCz = dCf(:,3);
            app.ddCx = ddCf(:,1); app.ddCy = ddCf(:,2); app.ddCz = ddCf(:,3);

            app.qx = app.A \ (x - app.Cx);
            app.qy = app.A \ (y - app.Cy);
            app.qz = app.A \ (z - app.Cz);

            app.TOF = (app.iPeriod + app.fPeriod)/2;

            z0 = [app.qx; app.qy; app.qz; app.TOF];
            

            xtest = app.A *app.qx + app.Cx;
            ytest = app.A *app.qy + app.Cy;
            ztest = app.A *app.qz + app.Cz;


            e = eye(length(z0));
            e = e(1,:);


            problem = struct('objective',@transferGenCostFcn,'x0',z0 ...
                ,'Aineq',e,'bineq', app.Tmax, 'Aeq', [], 'beq', [] ...
                ,'lb', [], 'ub', [], 'nonlcon', [] ...
            ,'solver','fmincon', 'options', fmc_optsfcn);

            % with nonlinear constraint
            % problem = struct('objective',@transferGenCostFcn,'x0',z0 ...
            %     ,'Aineq',e,'bineq', app.Tmax, 'Aeq', [], 'beq', [] ...
            %     ,'lb', [], 'ub', [], 'nonlcon', @(z) confunc(z,0.01866,0.00478) ...
            % ,'solver','fmincon', 'options', []);
            [z_opt, J_opt] = fmincon(problem);
            
            z = z_opt;
            app.TOF = z(end);
            N = 2*app.nffs+1;
            app.qx = z(1:N);
            app.qy = z(N+1:2*N);
            app.qz = z(2*N+1: 3*N);

            x = app.A * app.qx + app.Cx;
            y = app.A * app.qy + app.Cy;
            z = app.A * app.qz + app.Cz;
            

            hold(app.UIAxes,'on')
            app.transferArc = plot3(app.UIAxes,x,y,z,'LineWidth',2, ...
                'Color',app.orange);
            app.FouriertransferArc = plot3(app.UIAxes, xtest,ytest,ztest ...
                ,'r--','LineWidth',2);
            hold(app.UIAxes,'off')
        end

    end

    % Cost function Definition
    function J = transferGenCostFcn(z)
        % z:  3(2nffs-3) x 1
        % = [TOF   qx'   qy'  qz']'
        % where qx = [a_x0 a_x3 b_x3 a_x4 b_x4 ... a_xN, b_xN]'
        %       qy = [a_y0 a_y3 b_y3 a_y4 b_y4 ... a_yN, b_yN]'
        %       qz = [a_z0 a_z3 b_z3 a_z4 b_z4 ... a_zN, b_zN]'
        
        app.TOF = z(end);
        N = 2*app.nffs+1;
        app.qx = z(1:N);
        app.qy = z(N+1:2*N);
        app.qz = z(2*N+1: 3*N);

        x = app.A * app.qx + app.Cx;
        y = app.A * app.qy + app.Cy;
        z = app.A * app.qz + app.Cz;

        xd = 1/app.TOF * (app.dA * app.qx + app.dCx);
        yd = 1/app.TOF * (app.dA * app.qy + app.dCy);
        zd = 1/app.TOF * (app.dA * app.qz + app.dCz);
        

        xdd = 1/app.TOF^2 * (app.ddA * app.qx + app.ddCx);
        ydd = 1/app.TOF^2 * (app.ddA * app.qy + app.ddCy);
        zdd = 1/app.TOF^2 * (app.ddA * app.qz + app.ddCz);

        for k = 1:app.m
            r1(k) = sqrt((x(k)-app.mu)^2 +y(k)^2 +z(k)^2);
            r2(k) = sqrt((x(k)-1+app.mu)^2 +y(k)^2 +z(k)^2);

            ux = xdd(k) - 2*yd(k) - x(k) +...
                (1-app.mu)*(x(k)+app.mu)/r1(k)^3 + ...
                app.mu*(x-1+app.mu)/r2(k)^3;

            uy = ydd(k) + 2*xd(k) - y(k) +...
                (1-app.mu)*y(k)/r1(k)^3 + app.mu*y(k)/r2(k)^3  ;

            
            uz = zdd(k) + (1-app.mu)*z(k)/r1(k)^3  + app.mu*z(k)/r2(k)^3;

            U = [ux;uy;uz];

            L(k) = U'*U;
        end

        J = trapz(t*app.TOF, L);
        sqrt(J)
    end


    function [c,ceq] = confunc(X, re, rm)
        c1 = re - r1;
        c2 = rm - r2;
        c = [c1;c2];
        ceq = [];
    end

end


