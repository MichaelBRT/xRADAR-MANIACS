function [state, idx, times] = getPoincare(sectionID, C, stability, orbitType)
% plotPoincare - Load and plot Poincare section data based on section ID and Jacobi constant.
%
% Inputs:
%    sectionID - A string indicating the Poincare section ('eX+', 'eY-', etc.)
%    C         - Jacobi constant (used as index into cell array)
%    stability - Boolean indicating if looking at stable or unstable direction
%    orbitType - Orbit family being plotted
%
% Outputs:
%    state - states of an orbit for the 2D Poincare section plot
%
% Path to .mat file
switch sectionID
    case 'eX'
        filePath = fullfile('Poincaré Section Data/eX/', [orbitType, '.mat']);
    case 'mX'
        filePath = fullfile('Poincaré Section Data/mX/', [orbitType, '.mat']);
    case 'eY'
        filePath = fullfile('Poincaré Section Data/eY/', [orbitType, '.mat']);
    case 'mY'
        filePath = fullfile('Poincaré Section Data/mY/', [orbitType, '.mat']);
end
if ~isfile(filePath)   ,    error('File not found: %s', filePath);    end
    
% Load data
S = load(filePath);
if stability  
    if ~isfield(S,'Ws_Section_Data') , error('Section_Data not found in %s.',filePath); end
    sectionData = S.Ws_Section_Data;
    sectionTime = S.ts_Section_Data;
else   
    if ~isfield(S,'Wu_Section_Data') , error('Section_Data not found in %s.',filePath); end
    sectionData = S.Wu_Section_Data;
    sectionTime = S.tu_Section_Data;
end
  
Cs = S.C;
[~, idx] = min(abs(Cs - C)); % determine index of relevant data via Jacobi constant
state = sectionData{idx};
times = sectionTime{idx};
    
end
