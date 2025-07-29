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

%  ADD TO THE RIGHT SIDE FO THE FIGURE (ANNOTATION COMMAND)
% FIRST BOX GIVE INFO ABOUT MANEUVER (DV, THRUST, TOF, INITIAL/FINAL ORBIT)
% INSTANTANEOUS INFO, (ON MANIFOLD, CURRENT THRUST, CURRENT CUMULATIVE DV, TIME, UX, UY, THETA)


% --------------------- FIGURE SETUP ------------------------

f = figure(); % Create figure window
movegui(f, 'center'); % Optional: centers window
ax = axes(f); % ensures ax is explicitly tied to the figure
% f.Theme = 'dark';
% ax = gca;     % Get current axes

% Set axis limits based on initial trajectory data
xlims = 1.1*[min(init_arc(:,1)),max(init_arc(:,1))];
ylims = 1.1*[min(init_arc(:,2)),max(init_arc(:,2))];

% Plot initial and final orbits
hold(ax,'on');
plot(ax,Xi(:,1),Xi(:,2),'LineWidth',2,'Color',blue)
plot(ax,Xf(:,1),Xf(:,2),'LineWidth',2,'Color',purple)
% Plot Earth and Moon system
utils.drawEarthMoonSystem(ax,1,Ci); 

xlim(xlims);
ylim(ylims);


% ----------------- ANIMATED OBJECTS SETUP ------------------

% Initialize handle for transfer trajectory line (orange)
h_traj = plot(ax,init_arc(1,1), init_arc(1,2), 'LineWidth', 1.5,'Color',orange);
% Initialize handle for moving spacecraft (white circle with orange edge)
h_obj = plot(ax, init_arc(1,1), init_arc(1,2), 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', 'w','MarkerEdgeColor',orange);
% Initialize handle for thrust vector (red arrow)
h_thrustVec = quiver(init_arc(1,1), init_arc(1,2),0,0,10 ...
    ,"filled",'LineWidth',3,'Color','r','MaxHeadSize',1);
if ~isvalid(h_traj) || ~isvalid(h_obj) || ~isvalid(h_thrustVec)
    error("One or more plot handles are invalid. Check initialization.");
end


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


% ---------------- STATIC ANNOTATION PANEL -------------------------
maneuverInfo = sprintf([ ...
    'Maneuver Info:\n' ...
    'Total ΔV: %.3f\n' ...
    'Max Thrust: %.3f\n' ...
    'Time of Flight: %.3f TU\n' ...
    'Initial Orbit: %s\n' ...
    'Target Orbit: %s'], ...
    deltaV_req, umax, tTU(end), initOrbitName, targetOrbitName);

annotation('textbox', [0.75 0.6 0.2 0.3], 'String', maneuverInfo, ...
    'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 10, 'FontWeight','bold');


% ---------------- DYNAMIC TEXT OBJECT (instantaneous info) ---------------
h_dynamicText = text(ax, xlims(1)+0.05*range(xlims), ylims(2)-0.05*range(ylims), '', ...
    'FontSize', 10, 'VerticalAlignment', 'top', 'FontWeight','bold', 'BackgroundColor','w');


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
    onManifold = u_interp(k) < 1e-6;  % adjust threshold as needed
    % Thrust angle (in degrees)
    theta = atan2d(uy_N(k), ux_N(k));  % handles quadrant properly

    % Instantaneous info text
    dynamicText = sprintf([ ...
        'Instantaneous Info:\n' ...
        'On Manifold: %s\n' ...
        'Thrust Mag: %.4f\n' ...
        'Time: %.2f TU\n' ...
        'θ (deg): %.1f'], ...
        string(onManifold), u_interp(k), t_interp(k), theta);

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


