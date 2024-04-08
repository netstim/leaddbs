function [settings,eleNum] = ea_get_oss_reco(options, settings)
% Get electrode reconstruction parameters in OSS-DBS format.
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    options     % Lead-DBS options for electrode reconstruction and stimulation
    settings    % parameters for OSS-DBS simulation
end

% Reload reco since we need to decide whether to use native or MNI coordinates.
coords_mm = ea_load_reconstruction(options);
settings.contactLocation = coords_mm;
eleNum = length(coords_mm); % Number of electrodes
conNum = options.elspec.numel; % Number of contacts per electrode

% Save both native and MNI space y and head markers for OSS-DBS
settings.yMarkerNative = nan(eleNum, 3);
settings.yMarkerMNI = nan(eleNum, 3);
settings.headNative = nan(eleNum, 3);
settings.headMNI = nan(eleNum, 3);
[markersNative, markersMNI] = ea_get_markers(options);
for i=1:eleNum
    if ~isempty(markersNative) && ~isempty(markersNative(i).y)
    	settings.yMarkerNative(i,:) = markersNative(i).y;
    end
    if ~isempty(markersMNI) && ~isempty(markersMNI(i).y)
    	settings.yMarkerMNI(i,:) = markersMNI(i).y;
    end
    if ~isempty(markersNative) && ~isempty(markersNative(i).head)
    	settings.headNative(i,:) = markersNative(i).head;
    end
    if ~isempty(markersMNI) && ~isempty(markersMNI(i).head)
    	settings.headMNI(i,:) = markersMNI(i).head;
    end
end

% Head
settings.Implantation_coordinate = nan(eleNum, 3);
for i=1:eleNum
    if ~isempty(coords_mm{i})
        settings.Implantation_coordinate(i,:) = coords_mm{i}(1,:);
    end
end

% Tail
settings.Second_coordinate = nan(eleNum, 3);
for i=1:eleNum
    if conNum == 1 % Exception for electrode with only one contact
        if options.native
            settings.Second_coordinate(i,:) = markersNative(i).tail;
        else
            settings.Second_coordinate(i,:) = markersMNI(i).tail;
        end
    elseif ~isempty(coords_mm{i})
        if contains(options.elmodel, 'DIXI D08')
            settings.Second_coordinate(i,:) = coords_mm{i}(4,:);
        else
            settings.Second_coordinate(i,:) = coords_mm{i}(end,:);
        end
    end
end

%% Helper function to get markers in bothe native and MNI space
function [markersNative, markersMNI] = ea_get_markers(options)
    options.native = 1;
    try
        [~, ~, markersNative] = ea_load_reconstruction(options);
    catch
        markersNative = [];
        fprintf('\n')
        warning('off', 'backtrace');
        warning('Failed to load native reconstruction!');
        warning('on', 'backtrace');
    end
    
    options.native = 0;
    try
        [~, ~, markersMNI] = ea_load_reconstruction(options);
    catch
        markersMNI = [];
        fprintf('\n')
        warning('off', 'backtrace');
        warning('Failed to load MNI reconstruction!');
        warning('on', 'backtrace');
    end