function saveTransferData(mu,app)
    tol = 1e-10;
    odeopts = odeset('RelTol',3*tol,'AbsTol',tol);
    dynamics = @(t,x) utils.pcr3bp(t,x,mu);

    Ci = app.initJacobi.Value;
    Cf = app.targetJacobi.Value;

    initOrbitFile = fullfile('filtered PlanarOrbitData/',app.initDropdown.Value);
    targetOrbitFile = fullfile('filtered PlanarOrbitData/',app.targetDropdown.Value);

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

    



end