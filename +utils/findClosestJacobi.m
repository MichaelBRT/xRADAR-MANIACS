function idx = findClosestJacobi(data, targetC)
    [~, idx] = min(abs(data(:,8) - targetC));
end