function animationApp()
    % Load data
    load('scenarios\transfer_data.mat')
    % Expects: U, Xi, Xf, init_arc, tTU, Ci, initOrbitName, targetOrbitName, deltaV_req

    % Colors
    blue   = [0.07, 0.62, 1.00];
    orange = [0.988, 0.38, 0];
    purple = [0.659, 0, 1];

    % Normalize thrust
    ux = U(:,1); uy = U(:,2); u_mag = U(:,3);
    umax = max(u_mag);
    ux_N = ux / umax;
    uy_N = uy / umax;

    % --- Interpolation for Animation ---
    N = 5000; num_frames = 500; step = N / num_frames;
    X = interp1(tTU, init_arc, linspace(0, tTU(end), N), "spline");
    ux_N = interp1(tTU, ux_N, linspace(0, tTU(end), N), "spline");
    uy_N = interp1(tTU, uy_N, linspace(0, tTU(end), N), "spline");
    t_interp = linspace(tTU(1), tTU(end), N);

    % --- UIFigure & Layout ---
    appFig = uifigure('Name','Orbit Transfer Animation','Position',[100 100 1200 600]);
    grid = uigridlayout(appFig, [3, 3]);   
    grid.ColumnWidth = {'1x', '3x', '2x'};
    grid.RowHeight = {'fit'};    % 0.7x
    grid.Padding = [0 0 0 0];
    grid.RowSpacing = 1;
    grid.ColumnSpacing = 10;


    % Left Panel (Text Info)
    leftPanel = uipanel(grid);
    leftPanel.Layout.Row = 2;
    leftPanel.Layout.Column = 1;
    leftLayout = uigridlayout(leftPanel, [2,1]);
    leftLayout.RowHeight = {'1x', '1x'};
    
    T_em = 2.361 * 10^6;  D = 3.850*10^5;  n = 2*pi/T_em;  mass = 1000; %kg
    nd_to_N = 1000 * mass * n^2 * D;
    maxThrustNewt = nd_to_N * umax;
    maneuverInfo = sprintf([ ...
        'Maneuver Info:\n' ...
        'Initial Orbit: %s\n' ...
        'Target Orbit: %s\n' ...
        'Time of Flight: %.3f TU\n' ...
        'Total ΔV: %.3f m/s\n' ...
        'Max Thrust: %.3f N\n'], ...
        initOrbitName, targetOrbitName, tTU(end), deltaV_req, maxThrustNewt);

    staticText = uitextarea(leftLayout, ...
        'Value', splitlines(maneuverInfo), ...
        'Editable', 'off', ...
        'FontSize', 12, 'FontWeight','bold');
    
    dynamicText = uitextarea(leftLayout, ...
        'Value', {''}, ...
        'Editable', 'off', ...
        'FontSize', 12, ...
        'FontWeight','bold');

    % Middle Panel (Animation)
    midPanel = uipanel(grid);
    midPanel.Layout.Row = 2;
    midPanel.Layout.Column = 2;

    % ax = uiaxes(midPanel);
    midLayout = uigridlayout(midPanel, [1,1]);
midLayout.RowHeight = {'1x'};
midLayout.ColumnWidth = {'1x'};
midLayout.Padding = [5 5 5 5];
midLayout.RowSpacing = 0;
midLayout.ColumnSpacing = 0;

ax = uiaxes(midLayout);
ax.Layout.Row = 1;
ax.Layout.Column = 1;

    hold(ax,'on');

    % Plot orbits
    xlims = 1.1 * [min(init_arc(:,1)), max(init_arc(:,1))];
    ylims = 1.1 * [min(init_arc(:,2)), max(init_arc(:,2))];
    plot(ax, Xi(:,1), Xi(:,2), 'LineWidth', 2, 'Color', blue);
    plot(ax, Xf(:,1), Xf(:,2), 'LineWidth', 2, 'Color', purple);
    utils.drawEarthMoonSystem(ax, 1, Ci);
    xlim(ax, xlims);
    ylim(ax, ylims);

    % Animated handles
    hold(ax,'on');
    h_traj = plot(ax, init_arc(1,1), init_arc(1,2), 'LineWidth', 1.5, 'Color', orange);
    h_obj = plot(ax, init_arc(1,1), init_arc(1,2), 'o', ...
        'MarkerSize', 5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', orange);
    h_thrustVec = quiver(ax, init_arc(1,1), init_arc(1,2), 0, 0, 10, ...
        'filled', 'LineWidth', 2, 'Color', 'r', 'MaxHeadSize', 1);

    % Right Panel (Thrust Profile)
    rightPanel = uipanel(grid);
    rightPanel.Layout.Row = 2;
    rightPanel.Layout.Column = 3;

    rightLayout = uigridlayout(rightPanel, [1,1]);
    rightLayout.RowHeight = {'1x'};       rightLayout.ColumnWidth = {'1x'};
    rightLayout.Padding = [5 5 5 5];      rightLayout.RowSpacing = 0;
    rightLayout.ColumnSpacing = 0;
    thrustAx = uiaxes(rightLayout);  
    thrustAx.Layout.Row = 1;        thrustAx.Layout.Column = 1;

    
    % Interpolated thrust curve
    thrustCurve = plot(thrustAx, t_interp, interp1(tTU, u_mag, t_interp, "spline"), ...
    'k', 'LineWidth', 1.5);  % black thrust curve
    % Moving thrust dot (black with cyan outline)
    h_thrustDot = plot(thrustAx, t_interp(1), u_mag(1), 'o', ...
    'MarkerSize', 6, ...
    'MarkerEdgeColor', 'c', ...
    'MarkerFaceColor', 'k');


    xlabel(thrustAx, 'Time [TU]', 'FontWeight','bold');
    ylabel(thrustAx, 'Thrust Magnitude' ,'FontWeight','bold');
    thrustAx.XGrid = 'on';
    thrustAx.YGrid = 'on'; 

    drawnow;  % ensure layout is computed
outerPos = appFig.OuterPosition;
innerExtent = grid.Position;
excessHeight = outerPos(4) - innerExtent(4);
appFig.Position(4) = innerExtent(4) + excessHeight + 100;  % adjust to content


    % --- Animation Loop ---
    for k = 1:step:N
        xk = X(1:k,1); yk = X(1:k,2);
        h_traj.XData = xk;
        h_traj.YData = yk;
        h_obj.XData = xk(end);
        h_obj.YData = yk(end);
        h_thrustVec.XData = xk(end);
        h_thrustVec.YData = yk(end);
        h_thrustVec.UData = ux_N(k);
        h_thrustVec.VData = uy_N(k);
        % Update thrust dot
        h_thrustDot.XData = t_interp(k);
        h_thrustDot.YData = interp1(tTU, u_mag, t_interp(k), "spline");


        % Update dynamic text
        u_interp = hypot(ux_N(k), uy_N(k));
        uinterpNewt = nd_to_N * u_interp;
        onManifold = u_interp < 1e-6;
        theta = onManifold * 0 + ~onManifold * atan2d(uy_N(k), ux_N(k));
        dynamicText.Value = sprintf([ ...
            'Instantaneous Info:\nTime: %.2f TU\nOn Manifold: %s\n' ...
            'Thrust Mag: %.4f N\nThrust Angle θ: %.1f deg'], ...
            t_interp(k), string(onManifold), uinterpNewt, theta);

        drawnow;
    end
end
