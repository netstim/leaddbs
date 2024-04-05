function [S,settings,activeSources] = ea_check_stimSources(options,S,settings)
% Check which sources are active and insert blank stimulation for unilateral cases.
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    options     % Lead-DBS options for electrode reconstruction and stimulation
    S           % Lead-DBS stimulation settings
    settings    % parameters for OSS-DBS simulation
end

conNum = options.elspec.numel; % number of contacts
coords_mm = ea_load_reconstruction(options);
eleNum = length(coords_mm); % number of electrodes

if eleNum == 1
    % if right electrode only, insert null Stim Protocol for the left
    S = ea_add_StimVector_to_S(S, zeros(1, conNum),1);
    eleNum = 2;  % will be always assumed as 2 downstream
elseif all(isnan(settings.Implantation_coordinate(1,:)))
    % check if the right electrode is missing (eleNum is already 2)
    S = ea_add_StimVector_to_S(S, zeros(1, conNum),0);
end

nActiveSources = zeros(2,1); % count sources for each side
activeSources = nan(2,4);
for side = 1:2
    for source_index = 1:4
        if S.amplitude{side}(source_index) ~= 0 && ~isnan(S.amplitude{side}(source_index))
            nActiveSources(side) = nActiveSources(side) + 1;
            activeSources(side,source_index) = source_index;
        end
    end

    if nActiveSources(side) > 1
        if settings.adaptive_threshold
            ea_warndlg('Adaptive Thresholding is not supported for MultiSource. Switching to machine.vatsettings.butenko_ethresh')
            settings.adaptive_threshold = 0;
        end
    end
end