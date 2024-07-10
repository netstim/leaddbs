function ea_create_target_dict(targetType,targetActivation,stimfolder,connectome)

% create a dictionary with simple target / avoidance tracts
if startsWith(connectome, 'Multi-Tract: ')
    connectome = strrep(connectome, 'Multi-Tract: ', '');
    % Get tract names
    TractConnFolder = [ea_getconnectomebase, 'dMRI_MultiTract'];
    connFolder = [TractConnFolder, filesep, connectome];
    tractNames = strrep(ea_regexpdir(connFolder, '\.mat$', 0), [connFolder, filesep], '');
    tractNames = strrep(tractNames, '.mat', '');
else
    tractNames = connectome;
end

for side = 1:2

    % if the connectome does not specify hemisphere, assume the tract needs
    % to be targeted / avoided by both leads
    if side == 1
        side_suffix = '_rh';
        side_name = '_right';
    else
        side_suffix = '_lh';
        side_name = '_left';
    end

    for t = 1:length(tractNames)

        tractName = tractNames{t};
        
        if contains(tractName,'_lh')  || contains(tractName,'_rh') || contains(tractName,'_left')  || contains(tractName,'_right')
            if (side == 2 && (contains(tractName,'_lh')  || contains(tractName,'_left'))) || (side == 1 && (contains(tractName,'_rh')  || contains(tractName,'_right')))
                % check if side and suffix match
    
                if strcmp(targetType(t),"SSE")
                    jsonDict.Soft_SE_dict.(genvarname(['avoidance', side_suffix])).(genvarname(tractName)) = [targetActivation(t), 1.0, 1];  % [target_rate, mean_val (not relevant), sum_val (can be used to emphasize an importance of the pathway)]
                elseif strcmp(targetType(t),"CSE") == 1
                    jsonDict.SE_dict.(genvarname(['complete_avoidance', side_suffix])).(genvarname(tractName)) = [targetActivation(t), 1.0, 1];
                else
                    jsonDict.profile_dict.(genvarname(['target', side_suffix])).(genvarname(tractName)) = [targetActivation(t), 1.0, 1];
                end
            else
                continue
            end
        else
            % optimization for both sides will use the tract
            if strcmp(targetType(t),"SSE")
                jsonDict.Soft_SE_dict.(genvarname(['avoidance', side_suffix])).(genvarname(tractName)) = [targetActivation(t), 1.0, 1];  % [target_rate, mean_val (not relevant), sum_val (can be used to emphasize an importance of the pathway)]
            elseif strcmp(targetType(t),"CSE") == 1
                jsonDict.SE_dict.(genvarname(['complete_avoidance', side_suffix])).(genvarname(tractName)) = [targetActivation(t), 1.0, 1];
            else
                jsonDict.profile_dict.(genvarname(['target', side_suffix])).(genvarname(tractName)) = [targetActivation(t), 1.0, 1];
            end   
        end
    end
end

% save to the same file
jsonText2 = jsonencode(jsonDict);
fid = fopen([stimfolder,filesep,'target_profiles.json'], 'w');
fprintf(fid, '%s', jsonText2)
fclose(fid);