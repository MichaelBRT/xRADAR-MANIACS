function plotSingleOrbitSurface(states, isVertical, app, color, label, stride, type)

    if isempty(states), return; end
    states = states(1:stride:end, :);
    plotType = type;

    if isVertical
        % y vs ydot (Ax1), y vs xdot (Ax2)
        plot(app.sectionAx1, states(:,2), states(:,4), plotType, 'Color', color, 'DisplayName', label);
        plot(app.sectionAx2, states(:,2), states(:,3), plotType, 'Color', color, 'DisplayName', label);
    else
        % x vs ydot (Ax1), x vs xdot (Ax2)
        plot(app.sectionAx1, states(:,1), states(:,4), plotType, 'Color', color, 'DisplayName', label);
        plot(app.sectionAx2, states(:,1), states(:,3), plotType, 'Color', color, 'DisplayName', label);
    end
end


%{
        3D
    if isempty(states), return; end
    if isVertical
        X = states(1:stride:end,2); Y = states(1:stride:end,3); Z = states(1:stride:end,4);
    else
        X = states(1:stride:end,1); Y = states(1:stride:end,3); Z = states(1:stride:end,4);
    end
    [xq, yq] = meshgrid(linspace(min(X), max(X), 30), linspace(min(Y), max(Y), 30));
    zq = griddata(X, Y, Z, xq, yq, 'natural');
    surf(app.poincareAx, xq, yq, zq, 'FaceAlpha', 0.9, 'EdgeColor', 'none', ...
         'FaceColor', color, 'DisplayName', label);
%}