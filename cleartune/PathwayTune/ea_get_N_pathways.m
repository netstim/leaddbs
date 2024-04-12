function N_tracts = ea_get_N_pathways(ActivationProfileDict)

    % get total number of pathways (tracts) exported from Fiber Filtering
    % ActivationProfileDict - json file generated when exporting the tracts

    jsonText = fileread(ActivationProfileDict);
    ProfileDict = jsondecode(jsonText); 

    N_tracts = 0;
    fn = fieldnames(ProfileDict.profile_dict);
    for count = 1:numel(fn)
        N_tracts = N_tracts + length(fieldnames(ProfileDict.profile_dict.(fn{count})));
    end

    if isfield(ProfileDict,'SE_dict')
        fn = fieldnames(ProfileDict.SE_dict);
        for count = 1:numel(fn)
            N_tracts = N_tracts + length(fieldnames(ProfileDict.SE_dict.(fn{count})));
        end
    end

    if isfield(ProfileDict,'Soft_SE_dict')
        fn = fieldnames(ProfileDict.Soft_SE_dict);
        for count = 1:numel(fn)
            N_tracts = N_tracts + length(fieldnames(ProfileDict.Soft_SE_dict.(fn{count})));
        end
    end
end