function app = PoincareMapNetwork()
% Creating transferGen programmatic app
clear; clc; close all;
%__________________________________________________________________________
% --- Parameters ---
params = struct('mu',[], 'C',[], 'DU',[], 'L_r',[], 'zvc', [],...
                'xL1',[], 'xL2',[], 'xL3',[], 'Re',[], 'Rm',[]);
params.mu = 0.01215058560962404;  mu  = params.mu;
params.C = 3.1;                   C   = params.C;
params.DU = 3.844e5;              DU  = params.DU; % km
params.L_r = 0.04;  % (DU) Lagrange point radius
params.xL1 = Lptpos(mu,1);     
params.xL2 = Lptpos(mu,2);   
params.xL3 = Lptpos(mu,3);  
params.Re = 6378 / DU;  
params.Rm = 1737 / DU;  
params.opts = odeset('RelTol',1e-10,'AbsTol',1e-11);
params.propfunc = @(t,x) utils.pcr3bp(t,x,mu);

isDragging = false;

% Things to save
init_arc =[];
U = [];
thrust = [];
deltaV_req = [];
td = [];
tTU = [];
thrust_arc = [];
XX0i = [];
XX0f = [];



%__________________________________________________________________________
% --- UI Setup ---
app.fig = uifigure('Name', 'Poincaré Map Network', ...
                   'Position', [100, 100, 1000, 600], ...
                   'WindowStyle', 'normal');
app.grid = uigridlayout(app.fig, [2, 3]);
app.grid.RowHeight = {'5x', 'fit'};
app.grid.ColumnWidth = {'3x', '2x', '1x'};

% Colors
app.blue     = [0.07, 0.62, 1.00];
app.orange   = [0.988, 0.38, 0];
app.purple   = [0.659, 0, 1];
app.gray     = [0.1, 0.1, 0.1];
app.cmap1    = [1.0000, 0.3756, 0.0000;1.0000, 0.4073, 0.2000;0.8000, 0.1302, 0.2401;0.8000, 0.5203, 0.0983;1.0000, 0.6638, 0.1206;1.0000, 0.0763, 0.0163;0.8000, 0.0476, 0.0423;0.8000, 0.3039, 0.0859;1.0000, 0.0000, 0.2211;1.0000, 0.2088, 0.3026];
app.cmap2    = [0.1048, 0.8849, 0.1224;0.4136, 0.1132, 0.8596;0.7403, 0.7580, 0.0848;0.1087, 0.9104, 0.4755;0.1471, 0.2529, 0.9282;0.9655, 0.9329, 0.1046;0.6917, 1.0000, 0.1043;0.4281, 0.8636, 0.0891;0.0827, 0.0346, 0.8865;0.8378, 0.1478, 0.9130];
app.fontsize = 16;
app.fontsize2= 8;

app.p = {};
%__________________________________________________________________________
% --- Main Plot Panel (Top Left) ---

% --- Grid setup ---
app.mainplotgrid = uigridlayout(app.grid,[2,1]);
app.mainplotgrid.RowHeight = {'10x','3x'};
app.mainplotgrid.Layout.Row = 1;
app.mainplotgrid.Layout.Column = 1;

% --- Synodic Frame Plot ---
app.ax = uiaxes(app.mainplotgrid);
app.ax.Layout.Row = 1;
app.ax.Layout.Column = 1;
title(app.ax, 'Cislunar Environment (Synodic Frame)', 'FontSize', app.fontsize);
xlabel(app.ax, '$x$ [DU]', 'Interpreter', 'latex', 'FontSize', app.fontsize);
ylabel(app.ax, '$y$ [DU]', 'Interpreter', 'latex', 'FontSize', app.fontsize);
xlim(app.ax, [-1.3 1.3]);
ylim(app.ax, [-1.1 1.1]);
axis(app.ax, 'equal');
grid(app.ax, "on");
hold(app.ax, 'on');

% --- Jacobi Constant Selection Figure
app.Cax = uiaxes(app.mainplotgrid);
app.Cax.Layout.Row = 2;
app.Cax.Layout.Column = 1;

app.Cax.PickableParts = 'none';
hold(app.Cax, 'on');
app.Ci_scatter = scatter(app.Cax,0,0,7.5,app.blue,'filled');
app.Cf_scatter = scatter(app.Cax,0,0,7.5,app.purple,'filled');
app.Ci_select_line = plot(app.Cax,[0,0],[0.5,2.5],'--','Color','k','LineWidth',1.5);
app.Ci_select_point = plot(app.Cax,0,2,'o','Color','k','MarkerFaceColor','k');

app.Cf_select_line = plot(app.Cax,[0,0],[0.5,2.5],'--','Color','k');
app.Cf_select_point = plot(app.Cax,0,1,'o','Color','k','MarkerFaceColor','k','LineWidth',1.5);
hold(app.Cax,'off');

%app.Ci_select_point

ylim(app.Cax, [0.5, 2+0.5]);
yticks(app.Cax, 1:2);
yticklabels(app.Cax, {'Target','Initial'});
xlabel(app.Cax, 'Jacobi Constant');
title(app.Cax, 'Drag, Click to Select C_J');


app.fig.WindowButtonDownFcn = @mouseDown;
app.fig.WindowButtonUpFcn = @mouseUp;



% --- Bottom Left Controls --- 
app.controlsRow = uipanel(app.grid);
app.controlsRow.Layout.Row = 2;
app.controlsRow.Layout.Column = 1;

controlsLayout = uigridlayout(app.controlsRow, [5, 4]);
controlsLayout.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit'};
controlsLayout.ColumnWidth = {'1x', '1x', '1x', 'fit'};


% --- Initial Orbit Controls ---
label1 = uilabel(controlsLayout, 'Text', 'Initial Orbit:');
label1.Layout.Row = 1; label1.Layout.Column = 1;

app.initDropdown = uidropdown(controlsLayout, 'Items', utils.getOrbitFileList());
app.initDropdown.Layout.Row = 1; app.initDropdown.Layout.Column = 2;

app.initJacobi = uieditfield(controlsLayout, 'numeric', 'Value', C);
app.initJacobi.Layout.Row = 1; app.initJacobi.Layout.Column = 3;

% --- Target Orbit Controls ---
label2 = uilabel(controlsLayout, 'Text', 'Target Orbit:');
label2.Layout.Row = 2; label2.Layout.Column = 1;

app.targetDropdown = uidropdown(controlsLayout, 'Items', utils.getOrbitFileList());
app.targetDropdown.Layout.Row = 2; app.targetDropdown.Layout.Column = 2;

app.targetJacobi = uieditfield(controlsLayout, 'numeric', 'Value', C);
app.targetJacobi.Layout.Row = 2; app.targetJacobi.Layout.Column = 3;

% --- Poincaré Section Controls ---
label3 = uilabel(controlsLayout, 'Text', 'Section:');
label3.Layout.Row = 3; label3.Layout.Column = 1;

app.sectionChoice = uidropdown(controlsLayout, ...
    'Items', {'Earth X', 'Moon X', 'Earth Y', 'Moon Y'}, ...
    'ItemsData', {'eX', 'mX', 'eY', 'mY'}, ...
    'Placeholder', 'Select section');
app.sectionChoice.Layout.Row = 3; app.sectionChoice.Layout.Column = 2;


%__________________________________________________________________________
% --- Poincaré Panel ---
app.poincarePanel = uipanel(app.grid, 'Title', 'Poincaré Section', 'FontSize', app.fontsize);
app.poincarePanel.Layout.Row = [1, 2];  % Span both rows
app.poincarePanel.Layout.Column = 2;
poincareLayout = uigridlayout(app.poincarePanel, [3,1]);
poincareLayout.RowHeight = {'4x', 'fit', '4x'};
poincareLayout.ColumnWidth = {'1x'};

app.sectionAx1 = uiaxes(poincareLayout);
app.sectionAx1.Layout.Row = 1;
title(app.sectionAx1, 'Section Plot 1', 'FontSize', app.fontsize);

app.sectionAx2 = uiaxes(poincareLayout);
app.sectionAx2.Layout.Row = 3;
title(app.sectionAx2, 'Section Plot 2', 'FontSize', app.fontsize);

% --- Slider Section ---
sliderRow = uigridlayout(poincareLayout, [1, 3]);
sliderRow.Layout.Row = 2;
sliderRow.Layout.Column = 1;
sliderRow.ColumnWidth = {25, 50, '1x'};
sliderRow.RowHeight = {'fit'};

tauLabel = uilabel(sliderRow, ...
    'Text', '$\tau$', ...
    'Interpreter', 'latex', ...
    'FontName', 'Segoe UI Symbol', ...
    'FontSize', app.fontsize + 12, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'center');  
tauLabel.Layout.Row = 1;
tauLabel.Layout.Column = 1;

% Slider
T_max = 10;  % Replace dynamically as needed
app.poincareSlider = uislider(sliderRow, ...
    'Limits', [0 T_max], ...
    'Value', 1, ...
    'MajorTicks', 0:2:T_max, ...
    'ValueChangingFcn', @(src, event) onSliderChanging(src, event), ...
    'ValueChangedFcn', @(src, event) onSliderReleased(src, event));
app.poincareSlider.Layout.Row = 1;
app.poincareSlider.Layout.Column = 3;
% Slider Val
app.sliderValueLabel = uilabel(sliderRow, ...
    'Text', sprintf('= %.2f', app.poincareSlider.Value), ...
    'FontName', 'Segoe UI Symbol', ...
    'FontSize', app.fontsize, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'center');
app.sliderValueLabel.Layout.Row = 1;
app.sliderValueLabel.Layout.Column = 2;


%__________________________________________________________________________
% --- Legend Panel ---

app.legendPanel = uipanel(app.grid, 'Title', 'Legends', 'FontSize', app.fontsize);
app.legendPanel.Layout.Row = [1, 2];
app.legendPanel.Layout.Column = 3;

legendLayout = uigridlayout(app.legendPanel, [3,1]);
legendLayout.RowHeight = {'1x', 2, '1x'};
legendLayout.ColumnWidth = {'1x'};

app.legendLayout1 = uigridlayout(legendLayout);
app.legendLayout1.Layout.Row = 1;
app.legendLayout1.Layout.Column = 1;
app.legendLayout1.RowHeight = {};
app.legendLayout1.ColumnWidth = {20, '1x'};

divider = uipanel(legendLayout, 'BackgroundColor', [0.6 0.6 0.6]);
divider.Layout.Row = 2;
divider.Layout.Column = 1;

app.legendLayout2 = uigridlayout(legendLayout);
app.legendLayout2.Layout.Row = 3;
app.legendLayout2.Layout.Column = 1;
app.legendLayout2.RowHeight = {};
app.legendLayout2.ColumnWidth = {20, '1x'};


%__________________________________________________________________________
% --- Action Buttons ---
% Plots Initial and Target Orbits
orbitBtn = uibutton(controlsLayout, 'Text', 'Plot Orbits', ...
    'ButtonPushedFcn', @(btn, event) applyOrbits( ...
        app.initDropdown.Value, app.initJacobi.Value, ...
        app.targetDropdown.Value, app.targetJacobi.Value, ...
        params, app));
orbitBtn.Layout.Row = [1 2];   orbitBtn.Layout.Column = 4;

% Plots Poincaré Sections
poincareBtn = uibutton(controlsLayout, 'Text', 'Plot Poincaré', ...
    'ButtonPushedFcn', @(btn, event) ...
        onPoincareApplyPreLoad(app.sectionChoice.Value, params, app));
poincareBtn.Layout.Row = 3;   poincareBtn.Layout.Column = 4;

% Resets the App
resetBtn = uibutton(controlsLayout, 'Text', 'Reset', ...
    'ButtonPushedFcn', @(btn, event) utils.resetPlot(app, params, C));
resetBtn.Layout.Row = 4;   resetBtn.Layout.Column = 4;

% Plots the Initial Guess Transfer Arc
guessBtn = uibutton(controlsLayout, 'Text', 'Guess Transfer', ...
    'ButtonPushedFcn', @(btn, event) plotGuess());
guessBtn.Layout.Row = 4;  % place it in a new row
guessBtn.Layout.Column = 1;  

% --- Save Scenario button
app.saveScenarioButton = uibutton(controlsLayout, ...
    "Text", "Save Scenario","ButtonPushedFcn",@saveTransferData);
app.saveScenarioButton.Layout.Row = 5;
app.saveScenarioButton.Layout.Column = 1;

% Switch
app.guessModeSwitch = uiswitch(controlsLayout, 'toggle');
app.guessModeSwitch.Items = {'ΔV', 'Thrust'};
app.guessModeSwitch.Value = 'ΔV';
app.guessModeSwitch.Layout.Row = 4;
app.guessModeSwitch.Layout.Column = 2;

% Animation App Button
app.animationButton = uibutton(controlsLayout,'Text','Visualize Transfer', ...
                              'ButtonPushedFcn',@animationAppGen);
app.animationButton.Layout.Row = 4;
app.animationButton.Layout.Column = 3;

%__________________________________________________________________________
% --- Plot Earth, Moon, Lagrange Points --- 
hold(app.ax, 'on');
utils.defaultCislunar(app, params);
% --- Plot ZVC --- 
utils.ZVC(C, mu, app, params);
% --- Plot Poincaré sections ---  
utils.defaultPoincare(app, params);

% --- Default Orbit Selections
init_file = 'Lyapunov (L1).csv';
target_file = 'Lyapunov (L2).csv';
init_C = 3.15;
target_C = 3.15;
applyOrbits(init_file,init_C,target_file,target_C,params,app);


%__________________________________________________________________________
%%     ~   { ~ FUNCTIONS ~ }   ~
%__________________________________________________________________________

function onPoincareApplyPreLoad(sectionID, params, app)

% --- Set up for Orbits Through Section ---

    % Struct for Poincaré Data
    app.poincareData = struct( ...
        'Initial', struct('S', [], 'T', [], 'Handles', []), ...
        'Target', struct('S', [], 'T', [], 'Handles', []), ...
        'Others', []);

    % Clear old Poincaré section line(s) from left panel
    delete(findall(app.ax, 'Tag', 'poincare'));
    % Re-plot only the selected section
    utils.selectedSection(app, params, sectionID);
    
    % Clear both axes on Poincaré panel
    pax1 = app.sectionAx1;    pax2 = app.sectionAx2;
    cla(pax1); cla(pax2);
    hold(pax1, 'on'); hold(pax2, 'on');
    grid(pax1, 'on'); grid(pax2, 'on');

    % Get numOrbits for Poincaré section
    numOrbits = sum(~[dir(fullfile(pwd, 'Poincaré Section Data', sectionID)).isdir]);
    % Determine section type
    isVertical = contains(sectionID, 'Y');
    % Create color map
    % cmap = utils.customColormap(numOrbits);
    % cmap2 = 1 - cmap + app;
    cmap1 = app.cmap1;
    cmap2 = app.cmap2;
    plotType = '.';



% --- Initial and Target Orbits ---

    p = gobjects(1,4);
    
    % Initial Orbit (Unstable Manifold)
    initialOrbit = erase(app.initDropdown.Value, '.csv');
    initialC = app.initJacobi.Value;
    [initialStates, initialTime, ~] = utils.getPoincare(sectionID, initialC, 0, initialOrbit);
    app.poincareData.Initial.S = initialStates;
    app.poincareData.Initial.T = initialTime;
    if isVertical
        p(1) = plot(pax1, initialStates(:,2), initialStates(:,4), plotType, 'Color', app.blue, 'DisplayName', [initialOrbit ' (unstable)']);
        p(2) = plot(pax2, initialStates(:,2), initialStates(:,3), plotType, 'Color', app.blue, 'DisplayName', [initialOrbit ' (unstable)']);
    else
        p(1) = plot(pax1, initialStates(:,1), initialStates(:,4), plotType, 'Color', app.blue, 'DisplayName', [initialOrbit ' (unstable)']);
        p(2) = plot(pax2, initialStates(:,1), initialStates(:,3), plotType, 'Color', app.blue, 'DisplayName', [initialOrbit ' (unstable)']);
    end
    app.poincareData.Initial.Handles = [p(1), p(2)];
    
    % Target Orbit (Stable Manifold)
    targetOrbit = erase(app.targetDropdown.Value, '.csv');
    targetC = app.targetJacobi.Value;
    [targetStates, targetTime, ~] = utils.getPoincare(sectionID, targetC, 1, targetOrbit);
    app.poincareData.Target.S = targetStates;
    app.poincareData.Target.T = targetTime;
    if isVertical
        p(3) = plot(pax1, targetStates(:,2), targetStates(:,4), plotType, 'Color', app.purple, 'DisplayName', [targetOrbit ' (stable)']);
        p(4) = plot(pax2, targetStates(:,2), targetStates(:,3), plotType, 'Color', app.purple, 'DisplayName', [targetOrbit ' (stable)']);
    else
        p(3) = plot(pax1, targetStates(:,1), targetStates(:,4), plotType, 'Color', app.purple, 'DisplayName', [targetOrbit ' (stable)']);
        p(4) = plot(pax2, targetStates(:,1), targetStates(:,3), plotType, 'Color', app.purple, 'DisplayName', [targetOrbit ' (stable)']);
    end
    app.poincareData.Target.Handles  = [p(3), p(4)];

    % Start with Visibilty of Initial and Target Turned off
    set(p(:), 'Visible', 'on');


    
    % --- Other Orbits Through Section ---

    files = dir(fullfile('Poincaré Section Data', sectionID, '*.mat'));
    [~, idx] = sort({files.name});
    files = files(idx);

    % Struct for Other Orbits
    app.poincareData.Others = struct([]);
    T_max = 0;
    plotRow = 2; % row 1 is initial/target, others start at row 2

    % Loop through Orbits in Section
    for o = 1:numOrbits  

        % Get current orbit
        [~, orbit] = fileparts(files(o).name);
        % Skip Initial and Target Orbits
        if strcmp(orbit, initialOrbit) || strcmp(orbit, targetOrbit) , continue , end

        % Get Data
        try
            [u_state, u_time, ~] = utils.getPoincare(sectionID, initialC, 0, orbit); % unstable state
            [s_state, s_time, ~] = utils.getPoincare(sectionID, initialC, 1, orbit); % stable state
        catch
            warning('Skipping orbit %s due to missing data: %s', orbit);  continue 
        end
        if size(u_state, 2) < 4,  continue,  end

        % Store data
        app.poincareData.Others(plotRow-1).Name = orbit;
        app.poincareData.Others(plotRow-1).Su  = u_state;
        app.poincareData.Others(plotRow-1).Ss  = s_state;
        app.poincareData.Others(plotRow-1).Tu  = u_time;
        app.poincareData.Others(plotRow-1).Ts  = s_time;

        % Track T_max
        T_max = max([T_max, u_time(end), s_time(end)]);

        % --- PLOTTING ---
        %    & store handles
        handles = gobjects(1,4);

        if isVertical
            handles(1) = plot(pax1, u_state(:,2), u_state(:,4), plotType, 'Color', cmap1(o,:), 'DisplayName', [orbit ' (unstable)']);
            handles(2) = plot(pax1, s_state(:,2), s_state(:,4), plotType, 'Color', cmap2(o,:), 'DisplayName', [orbit ' (stable)']);
            handles(3) = plot(pax2, u_state(:,2), u_state(:,3), plotType, 'Color', cmap1(o,:), 'DisplayName', [orbit ' (unstable)']);
            handles(4) = plot(pax2, s_state(:,2), s_state(:,3), plotType, 'Color', cmap2(o,:), 'DisplayName', [orbit ' (stable)']);
            set(handles(:), 'Visible', 'off');
            xlabel(pax1, '$y$', 'Interpreter', 'latex', 'FontSize', app.fontsize);
            ylabel(pax1, '$\dot{y}$', 'Interpreter', 'latex', 'FontSize', app.fontsize);
            xlabel(pax2, '$y$', 'Interpreter', 'latex', 'FontSize', app.fontsize);
            ylabel(pax2, '$\dot{x}$', 'Interpreter', 'latex', 'FontSize', app.fontsize);
        else
            handles(1) = plot(pax1, u_state(:,1), u_state(:,4), plotType, 'Color', cmap1(o,:), 'DisplayName', [orbit ' (unstable)']);
            handles(2) = plot(pax1, s_state(:,1), s_state(:,4), plotType, 'Color', cmap2(o,:), 'DisplayName', [orbit ' (stable)']);
            handles(3) = plot(pax2, u_state(:,1), u_state(:,3), plotType, 'Color', cmap1(o,:), 'DisplayName', [orbit ' (unstable)']);
            handles(4) = plot(pax2, s_state(:,1), s_state(:,3), plotType, 'Color', cmap2(o,:), 'DisplayName', [orbit ' (stable)']);
            set(handles(:), 'Visible', 'off');
            xlabel(pax1, '$x$', 'Interpreter', 'latex', 'FontSize', app.fontsize);
            ylabel(pax1, '$\dot{y}$', 'Interpreter', 'latex', 'FontSize', app.fontsize);
            xlabel(pax2, '$x$', 'Interpreter', 'latex', 'FontSize', app.fontsize);
            ylabel(pax2, '$\dot{x}$', 'Interpreter', 'latex', 'FontSize', app.fontsize);
        end
        % Store Handles for Other Orbits
        app.poincareData.Others(plotRow-1).Handles = handles;

    end
    % Check Valid T_max
    if isempty(T_max) || ~isfinite(T_max) || T_max <= 0 ,    T_max = 1; end
    app.poincareSlider.Limits = [0 T_max];
    app.poincareSlider.MajorTicks = linspace(0, T_max, 6);



% --- LEGEND ---
    % Clear previous legend content
    delete(app.legendLayout1.Children);
    delete(app.legendLayout2.Children);
    
    % Get plot handles and flip for visual stack order
    ax1Children = flipud(findobj(pax1, 'Type', 'Line'));
    ax2Children = flipud(findobj(pax2, 'Type', 'Line'));
    


    % Collect all DisplayNames for Ax1
    orbitNames1 = unique(cellfun(@(s) erase(s, {' (unstable)',' (stable)'}), ...
                        get(ax1Children, 'DisplayName'), 'UniformOutput', false));
    for i = 1:numel(orbitNames1)
        orbitName = orbitNames1{i};
        % Find all lines whose DisplayName contains this orbitName
        hLines = ax1Children(contains(get(ax1Children, 'DisplayName'), orbitName));
        % Create a colored panel from the first line
        q = uipanel(app.legendLayout1, 'BackgroundColor', [1 1 1], 'BorderType', 'none');
        q.Layout.Row = i;
        q.Layout.Column = 1;
        q.Layout.Row = i;
        % Split: left = unstable, right = stable
        g = uigridlayout(q, [1,2]);  g.Padding = [0 0 0 0];         g.RowSpacing = 0;  
        g.ColumnSpacing = 0;         g.ColumnWidth = {'1x', '1x'};  g.RowHeight = {app.fontsize2 + 4}; % small swatches
        % Find handles for unstable/stable
        u = hLines(contains(get(hLines, 'DisplayName'), '(unstable)'));
        s = hLines(contains(get(hLines, 'DisplayName'), '(stable)'));
        if ~isempty(u), uipanel(g,'BackgroundColor',u(1).Color);
        else, uipanel(g, 'BackgroundColor', app.gray);  end
        if ~isempty(s), uipanel(g, 'BackgroundColor', s(1).Color);
        else, uipanel(g, 'BackgroundColor', app.gray);  end

        % Default visible if any line is visible
        isVisible = any(arrayfun(@(h) strcmp(h.Visible, 'on'), hLines));
        cb = uicheckbox(app.legendLayout1, ...
            'Text', orbitName, ...
            'Value', isVisible, ...
            'FontSize', app.fontsize2, ...
            'ValueChangedFcn', @(cb, ~) set(hLines, 'Visible', logical(cb.Value)));
        cb.Layout.Row = i;
        cb.Layout.Column = 2;
    end
    app.legendLayout1.RowHeight = repmat({22}, 1, numel(orbitNames1));



    % Collect all DisplayNames for Ax2
    orbitNames2 = unique(cellfun(@(s) erase(s, {' (unstable)',' (stable)'}), ...
                        get(ax2Children, 'DisplayName'), 'UniformOutput', false));

    for i = 1:numel(orbitNames2)
        orbitName = orbitNames2{i};    
        hLines = ax2Children(contains(get(ax2Children, 'DisplayName'), orbitName));
        q = uipanel(app.legendLayout2, 'BackgroundColor', [1 1 1], 'BorderType', 'none');
        q.Layout.Row = i;
        q.Layout.Column = 1;
        q.Layout.Row = i;
        % Split: left = unstable, right = stable
        g = uigridlayout(q, [1,2]);  g.Padding = [0 0 0 0];         g.RowSpacing = 0;  
        g.ColumnSpacing = 0;         g.ColumnWidth = {'1x', '1x'};  g.RowHeight = {app.fontsize2 + 4};
        % Find handles for unstable/stable
        u = hLines(contains(get(hLines, 'DisplayName'), '(unstable)'));
        s = hLines(contains(get(hLines, 'DisplayName'), '(stable)'));
        if ~isempty(u), uipanel(g,'BackgroundColor',u(1).Color);
        else, uipanel(g, 'BackgroundColor', app.gray);  end
        if ~isempty(s), uipanel(g, 'BackgroundColor', s(1).Color);
        else, uipanel(g, 'BackgroundColor', app.gray);  end

        % Default visible if any line is visible
        isVisible = any(arrayfun(@(h) strcmp(h.Visible, 'on'), hLines));
        cb = uicheckbox(app.legendLayout2, ...
            'Text', orbitName, ...
            'Value', isVisible, ...
            'FontSize', app.fontsize2, ...
            'ValueChangedFcn', @(cb, ~) set(hLines, 'Visible', logical(cb.Value)));
        cb.Layout.Row = i;
        cb.Layout.Column = 2;
    end
    app.legendLayout2.RowHeight = repmat({22}, 1, numel(orbitNames2));
    
end


%__________________________________________________________________________

function applyOrbits(initFile, initC, targetFile, targetC, params, app)
    folder = 'filtered PlanarOrbitData';
    ax = app.ax;
    % Delete old orbits
    C = initC;
    utils.resetPlot(app, params, C);
    % delete(findall(ax, 'Tag', 'orbit'));
    
    hold(ax, 'on');

    % Plot ZVC
    utils.ZVC(initC, params.mu, app);
    
    % Load and plot initial orbit
    data1 = readmatrix(fullfile(folder, initFile));
    idx1 = utils.findClosestJacobi(data1, initC);
    x01 = data1(idx1, 2:3);
    v01 = data1(idx1, 5:6);
    C01 = data1(idx1,8);
    T1 = data1(idx1, 9);
    [~, S1] = ode113(@(t,x) CR3BPMC2D(x, mu), [0 T1], [x01 v01], ...
                     params.opts);
    h1 = plot(ax, S1(:,1), S1(:,2), 'Color', app.blue, 'LineWidth', 1.5);
    utils.updateC(C01,1,app);
    set(h1, 'Tag', 'orbit');    

    % Load and plot target orbit
    data2 = readmatrix(fullfile(folder, targetFile));
    idx2 = utils.findClosestJacobi(data2, targetC);
    x02 = data2(idx2, 2:3);
    v02 = data2(idx2, 5:6);
    C02 = data2(idx2,8);
    T2 = data2(idx2, 9);
    [~, S2] = ode113(@(t,x) CR3BPMC2D(x, mu), [0 T2], [x02 v02], ...
                     params.opts);
    h2 = plot(ax, S2(:,1), S2(:,2), 'Color', app.purple, 'LineWidth', 1.5);
    utils.updateC(C02,2,app);
    set(h2, 'Tag', 'orbit');


    % Plot Jacobi Constant Range and Cursors
    C1 = data1(:,8);
    C2 = data2(:,8);

    app.Ci_scatter.XData = C1';
    app.Ci_scatter.YData = 2*ones(size(C1));

    app.Cf_scatter.XData = C2';
    app.Cf_scatter.YData = ones(size(C2));

end

%__________________________________________________________________________
function plotGuess()
    % Theoretical min DV calc:
    dVmin = utils.theoretical_min_dV(app.initJacobi.Value,app.targetJacobi.Value,mu);
    
    % --- Access handles
    ax = app.ax;
    pax1 = app.sectionAx1;
    pax2 = app.sectionAx2;

    % --- 
    [init_arc, U, thrust, deltaV_req,td,tTU, thrust_arc,XX0i,XX0f] = utils.get_initial_arc(mu,app);
    

    if ~isempty(init_arc)
    % --- Plot guess trajectory in main synodic plot
        hold(ax,'on')
        plot(ax, init_arc(:,1),init_arc(:,2), '--', 'Color', app.orange, ...
            'LineWidth', 1.5, 'DisplayName', 'Guess Trajectory');
        plot(ax,thrust_arc(:,1),thrust_arc(:,2),'LineWidth',1.5, ...
            'Color','r');
        plot(ax,[thrust_arc(1,1), thrust_arc(end,1)], ...
            [thrust_arc(1,2),thrust_arc(end,2)],'r','LineStyle','none' ...
            ,'MarkerFaceColor','w','LineWidth',1.5,'Marker','o');
        fprintf('delta V = %f \n max Thrust = %f \n ',deltaV_req,max(thrust))
        hold(ax,'off')


% --- Plot Thrust Profile

        % figure()
        % plot(td,thrust,'LineWidth',2)
        % thr_ax = gca;
        % grid('on')
        % thr_ax.Units = 'normalized';
        % % Get full figure-space position of plot box (not just axes container)
        % annotation_text = {['$\Delta V = $', num2str(1000*deltaV_req) ' m/s'],...
        % ['$\max \mathcal{T} = $' num2str(max(thrust)), ' N'],...
        % ['Time of Flight: ', num2str(max(td)),' days'], ...
        % ['Theoretical Minimum $\Delta V$: ',num2str(dVmin),' m/s']};
        %  axis tight
        % xl = xlim;
        % yl = ylim;
        % text(xl(1), yl(2), annotation_text, ...
        %      'HorizontalAlignment', 'left', ...
        %      'VerticalAlignment', 'top', ...
        %      'FontSize', 12, ...
        %      'Interpreter','latex', ...
        %      'Color','k')
        % xlabel('Time $t$ [days]','FontSize',14,'Interpreter','latex')
        % ylabel('Thrust $\mathcal{T}$ [N]','FontSize',14,'Interpreter','latex')
        % title('Thrust Profile of Transfer','Interpreter','latex')
       
    % --- Plot initial & target points on Poincaré plots (example values)


        hold(pax1,"on")
        plot(pax1, XX0i(2), XX0i(4), 'o', ...
            'MarkerFaceColor', app.blue, 'MarkerEdgeColor', 'k', ...
            'DisplayName', 'Initial Point');
    
        plot(pax1, XX0f(2), XX0f(4), 'o', ...
            'MarkerFaceColor', app.purple, 'MarkerEdgeColor', 'k', ...
            'DisplayName', 'Target Point');
        hold(pax1,"off")
    
        hold(pax2,"on")
        plot(pax2, XX0i(2), XX0i(3), 'o', ...
            'MarkerFaceColor', app.blue, 'MarkerEdgeColor', 'k', ...
            'DisplayName', 'Initial Point');
    
        plot(pax2, XX0f(2), XX0f(3), 'o', ...
            'MarkerFaceColor', app.purple, 'MarkerEdgeColor', 'k', ...
            'DisplayName', 'Target Point');
        hold(pax2,"off")
    end
end

% Jacobi Constant Slider
%______________________________________________________________________
   function mouseDown(~, ~)
        cp = app.Cax.CurrentPoint;
        xClick = cp(1,1);
        yClick = cp(1,2);

        % Extract current Jacobi constants
        Ci = app.Ci_select_point.XData;
        Cf = app.Cf_select_point.XData;


        di = sqrt((xClick - Ci)^2 + (yClick - 2)^2);
        df = sqrt((xClick - Cf)^2 + (yClick - 1)^2);

        if di < df
            index = 1;
        else
            index = 2;
        end

        C_vals = [Ci, Cf];
        y_vals = [2,1];

        % Decide which value
        if abs(xClick - C_vals(index)) < 0.01 && yClick >= 0 && yClick <= 2.5
            isDragging = true;
            app.fig.WindowButtonMotionFcn = @(src,event) mouseMove(src,event,index);
        end
    end

    % --- Mouse moved while dragging ---
    function mouseMove(~, ~,index)
        if isDragging
            C_new = app.Cax.CurrentPoint(1,1);
            utils.updateC(C_new,index,app);
        end
    end

    % --- Mouse released ---
    function mouseUp(~, ~)
        isDragging = false;
        app.fig.WindowButtonMotionFcn = '';
    end



%__________________________________________________________________________
function onSliderChanging(src, event)
    % Extract the live value as the slider is moving
    value = event.Value;
    app.sliderValueLabel.Text = sprintf('= %.2f', value);
    updateTrimmedPlots(value);
end

function onSliderReleased(src, ~)
    t = src.Value;
    app.sliderValueLabel.Text = sprintf('= %.2f', t);
    updateTrimmedPlots(t);
end


%__________________________________________________________________________
function updateTrimmedPlots(t)
    % t: current slider value

    % Exit early if data isn't loaded yet
    if ~isfield(app, 'poincareData') || isempty(app.poincareData), return  
    end
    isVertical = contains(app.sectionChoice.Value, 'Y');  % or however you store it

    % Trim Initial
    if ~isempty(app.poincareData.Initial.S)
        inds = find(app.poincareData.Initial.T <= t);
        trimInitial = app.poincareData.Initial.S(inds, :);
        if isVertical
            set(app.poincareData.Initial.Handles(1), 'XData', trimInitial(:,2), 'YData', trimInitial(:,4)); % Ax1
            set(app.poincareData.Initial.Handles(2), 'XData', trimInitial(:,2), 'YData', trimInitial(:,3)); % Ax2
        else
            set(app.poincareData.Initial.Handles(1), 'XData', trimInitial(:,1), 'YData', trimInitial(:,4)); % Ax1
            set(app.poincareData.Initial.Handles(2), 'XData', trimInitial(:,1), 'YData', trimInitial(:,3)); % Ax2
        end
    end

    % Trim Target
    if ~isempty(app.poincareData.Target.S)
        inds = find(app.poincareData.Target.T <= t);
        trimTarget = app.poincareData.Target.S(inds, :);
        if isVertical
            set(app.poincareData.Target.Handles(1), 'XData', trimTarget(:,2), 'YData', trimTarget(:,4)); % Ax1
            set(app.poincareData.Target.Handles(2), 'XData', trimTarget(:,2), 'YData', trimTarget(:,3)); % Ax2
        else
            set(app.poincareData.Target.Handles(1), 'XData', trimTarget(:,1), 'YData', trimTarget(:,4)); % Ax1
            set(app.poincareData.Target.Handles(2), 'XData', trimTarget(:,1), 'YData', trimTarget(:,3)); % Ax2
        end
    end
    disp(['Target has ' num2str(length(inds)) ' points'])

    % Loop over Other orbits
    for o = 1:numel(app.poincareData.Others)
        % Unstable
        indsU = find(app.poincareData.Others(o).Tu <= t);
        trimU = app.poincareData.Others(o).Su(indsU, :);
        if isVertical
            set(app.poincareData.Others(o).Handles(1), 'XData', trimU(:,2), 'YData', trimU(:,4)); % Ax1
            set(app.poincareData.Others(o).Handles(3), 'XData', trimU(:,2), 'YData', trimU(:,3)); % Ax2
        else
            set(app.poincareData.Others(o).Handles(1), 'XData', trimU(:,1), 'YData', trimU(:,4)); % Ax1
            set(app.poincareData.Others(o).Handles(3), 'XData', trimU(:,1), 'YData', trimU(:,3)); % Ax2
        end
        % Stable
        indsS = find(app.poincareData.Others(o).Ts <= t);
        trimS = app.poincareData.Others(o).Ss(indsS, :);
        if isVertical
            set(app.poincareData.Others(o).Handles(2), 'XData', trimS(:,2), 'YData', trimS(:,4)); % Ax1
            set(app.poincareData.Others(o).Handles(4), 'XData', trimS(:,2), 'YData', trimS(:,3)); % Ax2
        else
            set(app.poincareData.Others(o).Handles(2), 'XData', trimS(:,1), 'YData', trimS(:,4)); % Ax1
            set(app.poincareData.Others(o).Handles(4), 'XData', trimS(:,1), 'YData', trimS(:,3)); % Ax2
        end
        
    end
end


%__________________________________________________________________________

function saveTransferData(~,~)
    [outputfile,location] = uiputfile('*.mat','Save transfer data','transfer_data.mat');
    outputfile_full = fullfile(location,outputfile);
   
    % Check for cancel
    if isequal(outputfile, 0)
        disp('User canceled the save.');
        return;
    end



    tol = 1e-10;
    odeopts = odeset('RelTol',3*tol,'AbsTol',tol);
    dynamics = @(t,x) utils.pcr3bp(t,x,mu);

    Ci = app.initJacobi.Value;
    Cf = app.targetJacobi.Value;

    initOrbitFile = fullfile('filtered PlanarOrbitData/',app.initDropdown.Value);
    targetOrbitFile = fullfile('filtered PlanarOrbitData/',app.targetDropdown.Value);

    [~,initOrbitName] = fileparts(initOrbitFile);
    [~,targetOrbitName] = fileparts(targetOrbitFile);

    initOrbitData = readmatrix(initOrbitFile);
    targetOrbitData = readmatrix(targetOrbitFile);

    idx_i = utils.findClosestJacobi(initOrbitData,Ci);
    idx_f = utils.findClosestJacobi(targetOrbitData,Cf);

    Ci = initOrbitData(idx_i,8);
    Cf = targetOrbitData(idx_f,8);

    Xi_0 = initOrbitData(idx_i, [2,3,5,6])';
    Xf_0 = targetOrbitData(idx_f, [2,3,5,6])';

    Ti = initOrbitData(idx_i,9);
    Tf = targetOrbitData(idx_f,9);

    [ti,Xi] = ode45(dynamics,[0,Ti],Xi_0,odeopts);
    [tf,Xf] = ode45(dynamics,[0,Tf],Xf_0,odeopts);

    section_ID = app.sectionChoice.Value;

    [Wu_init,~,tu_init] = utils.getPoincare(section_ID,Ci,0,initOrbitName);
    [Ws_target,ts_target] = utils.getPoincare(section_ID,Cf,1,targetOrbitName);


    save(outputfile_full,'Ci','Cf',"Xi_0",'Xf_0','XX0i','XX0f',"Wu_init","Ws_target", ...
        "tu_init","ts_target","Xf","Xi","tf","ti","init_arc","U","deltaV_req" ...
        ,"thrust_arc","thrust","td", "tTU","initOrbitName","targetOrbitName" ...
        ,"section_ID");



end


%__________________________________________________________________________

function animationAppGen(~,~)
    % Load data
    Xi=XX0i; Xf=XX0f;
    %load('scenarios\transfer_data.mat')
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
    % Time vector for animation
    N = 5000; num_frames = 500; step = N / num_frames;
    t_interp = linspace(min(tTU), max(tTU), N);
    
    % Interpolated states and control
    X = interp1(tTU, init_arc, t_interp, "spline", 'extrap');
    ux_N = interp1(tTU, ux / umax, t_interp, "spline", 'extrap');
    uy_N = interp1(tTU, uy / umax, t_interp, "spline", 'extrap');
    u_interp_mag = interp1(tTU, u_mag, t_interp, "spline", 'extrap');


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

    hold(thrustAx,'on');
    % Moving thrust dot (black with cyan outline)
    h_thrustDot = plot(thrustAx, t_interp(1), u_interp_mag(1), 'o', ...
    'MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','c');

    xlabel(thrustAx, 'Time [TU]', 'FontWeight','bold');
    ylabel(thrustAx, 'Thrust Magnitude' ,'FontWeight','bold');
    thrustAx.XGrid = 'on';
    thrustAx.YGrid = 'on'; 

    drawnow;  % ensure layout is computed
    outerPos = appFig.OuterPosition;
    innerExtent = grid.Position;
    excessHeight = outerPos(4) - innerExtent(4);
    appFig.Position(4) = innerExtent(4) + excessHeight + 100;  % adjust to content

    % Interpolated thrust curve
    thrustCurve = plot(thrustAx, t_interp, u_interp_mag, ...
    'k', 'LineWidth', 1.5);
    % Set axis limits so it doesn't auto-zero
    xlim(thrustAx, [min(t_interp), max(t_interp)]);
    ylim(thrustAx, [0, 1.1 * max(u_interp_mag)]);


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
        h_thrustDot.YData = u_interp_mag(k);



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







end
