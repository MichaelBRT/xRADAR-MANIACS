clear, clc, close all;
load('scenarios\transfer_data.mat')
% U (control), Xi (initial orbit), Xf (final orbit), 
% init_arc (guess trajectory), tTU (time), Ci (initial Jacobi)


blue     = [0.07, 0.62, 1.00];
orange   = [0.988, 0.38, 0];
purple   = [0.659, 0, 1];
gray     = [0.1, 0.1, 0.1];


% Set scaling factor for visualizing thrust vector
control_scale_factor = 1e-2;
% Extract control input components from matrix U
ux = U(:,1);    % x-direction thrust
uy = U(:,2);    % y-direction thrust
u_mag = U(:,3); % magnitude of thrust
% Normalize control inputs for plotting
umax = max(u_mag);
ux_N = ux/umax;
uy_N = uy/umax;


% Set up video writer to save animation
savePath = fullfile('orbit animations/','orbit_animation_test.mp4');
v = VideoWriter(savePath, 'MPEG-4');
v.Quality = 100;
v.FrameRate = 30;
% Open the video file for writing
open(v); 


% --------------------- FIGURE SETUP ------------------------

f = figure('Position', [100, 100, 1000, 600]); 
t = tiledlayout(f, 1, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---- Left Panel (Orbit Plot) ----
panel_plot = uipanel(f, 'Title', '', 'BorderType', 'line');

ax = axes(panel_plot);
set(ax, 'Position', [0.1, 0.1, 0.85, 0.85]);
pbaspect(ax, [4 3 1]);
set(f, 'Resize', 'off'); 

% Set axis limits based on initial trajectory dat
x = init_arc(:,1);                    y = init_arc(:,2);
x_range = max(x) - min(x);            y_range = max(y) - min(y);
x_pad = 0.1 * x_range;                y_pad = 0.1 * y_range;
xlims = [min(x)-x_pad, max(x)+x_pad]; ylims = [min(y)-y_pad, max(y)+y_pad];
xlim(ax, xlims);                      ylim(ax, ylims);

hold(ax, 'on')
% Plot initial and final orbits
plot(ax,Xi(:,1),Xi(:,2),'LineWidth',2,'Color',blue)
plot(ax,Xf(:,1),Xf(:,2),'LineWidth',2,'Color',purple)
% Plot Earth and Moon system
utils.drawEarthMoonSystem(ax,1,Ci); 


% ----------------- ANIMATED OBJECTS SETUP ------------------

hold(ax,'on');
% Initialize handle for transfer trajectory line (orange)
h_traj = plot(ax,init_arc(1,1), init_arc(1,2), 'LineWidth', 1.5,'Color',orange);
% Initialize handle for moving spacecraft (white circle with orange edge)
h_obj = plot(ax, init_arc(1,1), init_arc(1,2), 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', 'w','MarkerEdgeColor',orange);
% Initialize handle for thrust vector (red arrow)
h_thrustVec = quiver(init_arc(1,1), init_arc(1,2),0,0,10 ...
    ,"filled",'LineWidth',3,'Color','r','MaxHeadSize',1);

% ---- Right Panel (Annotations) ----
panel_text = uipanel(f, 'Title', '', 'BorderType', 'line');
panel_text.Position = [0.76 0.1 0.22 0.85];

ax_ann = axes(panel_text);
axis(ax_ann, 'off')

% ------------------- INTERPOLATION --------------------------

N = 5000;            % Total interpolation steps
num_frames = 500;    % Desired number of video frames
% Step size in index for each frame
step = N/num_frames; % keep <500 frames
% Interpolate position data along the trajectory
X = interp1(tTU,init_arc,linspace(0,tTU(end),N),"spline");
% Interpolate normalized control inputs along the same time vector
ux_N = interp1(tTU,ux_N,linspace(0,tTU(end),N),"spline");
uy_N = interp1(tTU,uy_N,linspace(0,tTU(end),N),"spline");
% Interpolate time along trajectory
t_interp = linspace(tTU(1),tTU(end),N);

% ---------------- STATIC ANNOTATION PANEL -------------------------
maneuverInfo = sprintf([ ...
    'Maneuver Info:\n' ...   
    'Initial Orbit: %s\n' ...
    'Target Orbit: %s\n' ...
    'Time of Flight: %.3f TU\n' ...
    'Total ΔV: %.3f\n' ...
    'Max Thrust: %.3f\n'], ...
    initOrbitName, targetOrbitName, tTU(end), deltaV_req, umax);

set(ax_ann, 'Position', [0 0 1 1])  % Fill panel fully

% Instantaneous (top) annotation
h_dynamicText = annotation(panel_text, 'textbox', [0.05 0.55 0.9 0.4], ...
    'String', '', 'FontSize', 10, 'FontWeight', 'bold', ...
    'EdgeColor', 'none', 'Interpreter', 'none', 'HorizontalAlignment', 'left');

% Static (bottom) annotation
h_staticText = annotation(panel_text, 'textbox', [0.05 0.05 0.9 0.4], ...
    'String', maneuverInfo, 'FontSize', 10, 'FontWeight', 'bold', ...
    'EdgeColor', 'none', 'Interpreter', 'none', 'HorizontalAlignment', 'left');


% -------------------- ANIMATION LOOP ------------------------

for k = 1:step:N

    % Extract current trajectory segment
    xk = X(1:k,1);
    yk = X(1:k,2);

    % Update trajectory line
    h_traj.XData = xk;
    h_traj.YData = yk;

    % Update spacecraft position
    h_obj.XData = xk(end);
    h_obj.YData = yk(end);

    % Update thrust vector position and direction
    h_thrustVec.XData = xk(end);
    h_thrustVec.YData = yk(end);
    h_thrustVec.UData = ux_N(k);
    h_thrustVec.VData = uy_N(k);

    % Determine if on manifold (thrust ≈ 0)
    u_interp = sqrt(h_thrustVec.UData^2 + h_thrustVec.VData^2);
    onManifold = u_interp < 1e-6;  % adjust threshold as needed
    % Thrust angle (in degrees)
    if onManifold, theta=0;
    else, theta = atan2d(uy_N(k), ux_N(k));
    end
    
    % Instantaneous info text
    dynamicText = sprintf([ ...
        'Instantaneous Info:\n' ...
        'Time: %.2f TU\n' ...
        'On Manifold: %s\n' ...
        'Thrust Mag: %.4f\n' ...
        'θ (deg): %.1f'], ...
        t_interp(k), string(onManifold), u_interp, theta);

    h_dynamicText.String = dynamicText;

    % Render frame
    drawnow;
    % Capture current frame and write to video
    frame = getframe(f);
    writeVideo(v, frame);
end

% ------------------- CLEANUP -------------------------------

hold(ax,'off');
close(v);   % Finalize and close video file
close(f);   % Close figure


