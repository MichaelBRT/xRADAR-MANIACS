function draw_shaded_circle(ax, center, radius, color, fillCircle)
    if nargin < 4, color = [0.3, 0.6, 1]; fillCircle = true; end
    theta = linspace(0, 2*pi, 50);
    xCirc = center(1) + radius * cos(theta);
    yCirc = center(2) + radius * sin(theta);
    if fillCircle
        fill(ax, xCirc, yCirc, color, 'EdgeColor', 'none');
    else
        plot(ax, xCirc, yCirc, 'Color', color, 'LineStyle', '-.', 'LineWidth', 0.5);
    end
end