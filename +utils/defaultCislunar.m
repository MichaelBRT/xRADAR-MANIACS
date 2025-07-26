function defaultCislunar(app, params)
    ax = app.ax;          mu  = params.mu;  
    Re  = params.Re;      Rm  = params.Rm;     L_r = params.L_r; 
    xL1 = params.xL1;     xL2 = params.xL2;    xL3 = params.xL3;
    
    utils.draw_shaded_circle(ax, [-mu, 0], Re, app.blue, true);   % Earth
    utils.draw_shaded_circle(ax, [1 - mu, 0], Rm, app.gray, true);   % Moon
    utils.draw_shaded_circle(ax, [xL1, 0], L_r, [0 0 0], false);           % L1
    utils.draw_shaded_circle(ax, [xL2, 0], L_r, [0 0 0], false);           % L2
    utils.draw_shaded_circle(ax, [xL3, 0], L_r, [0 0 0], false);           % L3
    utils.draw_shaded_circle(ax, [0.5-mu,  sqrt(3)/2], L_r, [0 0 0], false);  % L4
    utils.draw_shaded_circle(ax, [0.5-mu, -sqrt(3)/2], L_r, [0 0 0], false);  % L5

end