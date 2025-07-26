function defaultPoincare(app, params)
    ax = app.ax;          mu  = params.mu;  
    Re  = params.Re;      Rm  = params.Rm;     L_r = params.L_r; 
    xL1 = params.xL1;     xL2 = params.xL2;          


    % Earth X
    plot(ax, [-1-5*L_r , -mu-Re*3], [0, 0], 'r--', 'LineWidth', 0.5, 'Tag', 'poincare'); 
    plot(ax, [-mu+Re*3, xL1], [0, 0], 'r--', 'LineWidth', 0.5, 'Tag', 'poincare');
    % Moon X
    plot(ax, [xL1, 1-mu-Rm*5], [0, 0], 'r--', 'LineWidth', 0.5, 'Tag', 'poincare');
    plot(ax, [1-mu+Rm*5, xL2+5*L_r], [0, 0], 'r--', 'LineWidth', 0.5, 'Tag', 'poincare');
    % Earth Y
    plot(ax, [-mu, -mu], [Re*3, 1], 'r--', 'LineWidth', 0.5, 'Tag', 'poincare'); 
    plot(ax, [-mu, -mu], [-1, -Re*3], 'r--', 'LineWidth', 0.5, 'Tag', 'poincare');
    % Moon Y
    plot(ax, [1-mu, 1-mu], [Rm*5, 0.5], 'r--', 'LineWidth', 0.5, 'Tag', 'poincare');
    plot(ax, [1-mu, 1-mu], [-0.5 , -Rm*5], 'r--', 'LineWidth', 0.5, 'Tag', 'poincare');
    

end