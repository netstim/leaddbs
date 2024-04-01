function [settings,runStatusMultiSource,activeSources] = ea_check_stimSources(S,settings,eleNum)

% Set stimSetMode flag
settings.stimSetMode = options.stimSetMode;
if settings.stimSetMode
    ea_warndlg("Not yet supported in V2")
    return
end

% if right electrode only, insert null Stim Protocol for the left
% this work around is not needed for the left only, handled by Lead-DBS
if eleNum == 1
    S = ea_add_StimVector_to_S(S, zeros(1, conNum),1);
    eleNum = 2;
elseif isempty(coords_mm)
    S = ea_add_StimVector_to_S(S, zeros(1, conNum),0);
end

% Initialize current control flag
settings.current_control = nan(eleNum, 1);

% Initalize Phi vector
settings.Phi_vector = nan(eleNum, conNum);

% Initialize grounding
settings.Case_grounding = zeros(eleNum, 1);

% Get the stimulation parameters from S in case stimSetMode is 0, otherwise
% they will be loaded directly from the Current_protocols_[0|1].csv files
% also get the center of the grid
settings.stim_center = nan(2, 3);

nActiveSources = zeros(2,1); % count sources for each side
activeSources = nan(2,4);
multiSourceMode = [0;0];
runStatusMultiSource = zeros(4,2);  % check status for each source
for side = 1:2
    for source_index = 1:4
        if S.amplitude{side}(source_index) ~= 0 && ~isnan(S.amplitude{side}(source_index))
            nActiveSources(side) = nActiveSources(side) + 1;
            activeSources(side,source_index) = source_index;
        end
    end

    if nActiveSources(side) > 1
        multiSourceMode(side) = 1;
        if settings.adaptive_threshold
            ea_warndlg('Adaptive Thresholding is not supported for MultiSource. Switching to machine.vatsettings.butenko_ethresh')
            settings.adaptive_threshold = 0;
        end
    end
end