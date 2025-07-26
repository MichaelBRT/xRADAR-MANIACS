function cmap = customColormap(N)
%CUSTOMCOLORMAP Generates N distinct RGB colors avoiding app.blue & app.purple
% Usage: cmap = utils.customColormap(50)
    
    if nargin < 1
        N = 50;
    end
    
    % Forbidden colors
    forbidden = [0.07, 0.62, 1.00;  % app.blue
        0.659, 0, 1.0];   % app.purple
    
    % HSV sampling
    numCandidates = 5 * N;  % oversample
    hues = linspace(0, 1, numCandidates)';
    sat = 0.75;  % Higher saturation for more pop
    val = 0.95;
    
    % Generate HSV candidates and convert to RGB
    hsvColors = [hues, repmat([sat, val], numCandidates, 1)];
    candidateRGB = hsv2rgb(hsvColors);
    
    % Prune near forbidden colors
    keep = true(numCandidates,1);
    minDist = 0.35;
    for i = 1:size(forbidden,1)
        dists = vecnorm(candidateRGB - forbidden(i,:), 2, 2);
        keep = keep & (dists > minDist);
    end
    candidateRGB = candidateRGB(keep, :);
    
    % Enforce distinctness from each other
    chosen = [];
    while size(chosen,1) < N && ~isempty(candidateRGB)
        chosen(end+1,:) = candidateRGB(1,:); %#ok<AGROW>
        % Remove all candidates too close to this one
        dists = vecnorm(candidateRGB - chosen(end,:), 2, 2);
        candidateRGB(dists < 0.2, :) = [];
    end
    
    % Fallback if not enough distinct
    if size(chosen,1) < N
        base = chosen;
        copies = ceil(N / size(base, 1));
        cmap = [];

        for i = 0:(copies - 1)
            % Slight hue rotation for each copy to reduce visual repetition
            hsvBase = rgb2hsv(base);
            hsvShifted = hsvBase;
            hsvShifted(:,1) = mod(hsvShifted(:,1) + rand * i, 1);  % shift hue by 15% of circle
            rgbShifted = hsv2rgb(hsvShifted);
            cmap = [cmap; rgbShifted]; %#ok<AGROW>
        end

    else
        cmap = chosen;
    end
    
    cmap = cmap(1:N, :);
end
