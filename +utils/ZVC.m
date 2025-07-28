function ZVC(C,mu,app,params)
    ax = app.ax;
    x = linspace(-1.5,1.5,1000);
    y = linspace(-1.5,1.5,1000);
    [X,Y] = meshgrid(x,y);

    Omega = @(x,y) 1/2*(x.^2 + y.^2) + (1-mu)./sqrt((x+mu).^2 + y.^2)+...
                mu./sqrt((x-1+mu).^2 + y.^2);

    Z = Omega(X,Y);
    params.zvc = Z;    

    hold(ax, 'on');

    % ZVC contour boundary at C/2
    Clevel = C / 2;
    %{
    contourData = contourc(x, y, Z, [Clevel Clevel]);

    % Extract contour polygon(s)
    idx = 1;
    while idx < size(contourData, 2)
        numPoints = contourData(2, idx);
        xPoly = contourData(1, idx+1:idx+numPoints);
        yPoly = contourData(2, idx+1:idx+numPoints);

        % Fill each forbidden region
        patch(ax, xPoly, yPoly, app.gray, ...
              'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
              'DisplayName', 'ZVC Forbidden', ...
              'Tag', 'ZVC');
        idx = idx + numPoints + 1;
    end
    %}

    % Draw the ZVC line on top
    contour(ax, X, Y, Z, [Clevel Clevel], ...
                  'LineColor', app.gray, 'LineWidth', 0.35, ...
                  'Tag', 'ZVC');
end

