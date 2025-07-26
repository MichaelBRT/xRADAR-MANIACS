function resetPlot(app, params, C)
% Reset Cislunar Axes
    ax = app.ax;
    cla(ax);
    hold(ax, 'on');
    axis(ax, 'equal');
    xlim(ax, [-1.3 1.3]);
    ylim(ax, [-1.1 1.1]);
    grid(ax, "on");

    % Reset both Poincaré section axes
    for ax = [app.sectionAx1, app.sectionAx2]
        cla(ax);
        hold(ax, 'on');
        grid(ax, 'on');
        title(ax, '', 'FontSize', app.fontsize);  % Clear title
        xlabel(ax, '', 'FontSize', app.fontsize);
        ylabel(ax, '', 'FontSize', app.fontsize);
    end

    % Clear interactive legend checkboxes and color panels
    delete(app.legendLayout1.Children);
    delete(app.legendLayout2.Children);

% Default Cislunar Environment Plot
    utils.defaultCislunar(app, params);
% Default ZVC 
    utils.ZVC(C, params.mu, app, params);
% Default Poincaré sections
    utils.defaultPoincare(app, params);



%{
           3D

    % Reset Poincaré section axes
    pax = app.poincareAx;
    cla(pax);
    hold(pax, 'on');
    view(pax, 3);  % Ensure it's in 3D mode
    grid(pax, 'on');
    xlabel(pax, '$x$', 'Interpreter', 'latex', 'FontSize', app.fontsize);
    ylabel(pax, '$\dot{x}$', 'Interpreter', 'latex', 'FontSize', app.fontsize);
    zlabel(pax, '$\dot{y}$', 'Interpreter', 'latex', 'FontSize', app.fontsize);
    title(pax, 'Section Plot', 'FontSize', app.fontsize);
    %}

end