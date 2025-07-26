function states = extractSectionCrossings(S, sectionID, params)
    tol = 1e-3;
    x = S(:,1); y = S(:,2);
    switch sectionID
        case {'eX', 'mX'}
            mask = abs(y) < tol;
        otherwise
            w = utils.getVerticalX(sectionID, params);
            mask = abs(x - w) < tol;
    end
    states = S(mask, :);
end