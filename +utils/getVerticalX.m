function xVal = getVerticalX(sectionID, params)
    mu = params.mu;
    switch sectionID
        case 'eY', xVal = -mu;
        case 'mY', xVal = 1 - mu;
        otherwise, xVal = 0;
    end
end