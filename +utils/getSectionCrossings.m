function states = getSectionCrossings(state0, T, N, sectionID, params)
    [~, S] = ode113(@(t,x) CR3BPMC2D(x, params.mu), [0 T*N], state0, params.opts);
    tol = 1e-3;
    if contains(sectionID, 'Y')
        mask = abs(S(:,1) - utils.getVerticalX(sectionID, params)) < tol;
    else
        mask = abs(S(:,2)) < tol;
    end
    states = S(mask, :);
end