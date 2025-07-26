function [crossingData, fileNames] = getOrbitsThroughSection(sectionID, initJacobi, params, numPeriods)
    folder = 'PlanarOrbitData';
    files = dir(fullfile(folder, '*.csv'));
    crossingData = {};
    fileNames = {};
    for i = 1:length(files)
        data = readmatrix(fullfile(folder, files(i).name));
        if isempty(data) || size(data, 2) < 9, continue; end
        idx = utils.findClosestJacobi(data, initJacobi);
        state0 = [data(idx,2:3), data(idx,5:6)];
        T = data(idx, 9);
        [~, S] = ode113(@(t,x) CR3BPMC2D(x, params.mu), [0 numPeriods*T], state0, params.opts);
        states = utils.extractSectionCrossings(S, sectionID, params);
        if ~isempty(states)
            crossingData{end+1} = states; %#ok<AGROW>
            fileNames{end + 1} = files(i).name; %#ok<AGROW>
        end
    end
end