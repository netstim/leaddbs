function ea_switch2bids

LeadRoot = ea_getearoot;

if ~isfolder(fullfile(LeadRoot, 'templates', 'space', 'MNI152NLin2009bAsym'))
    disp('Copy "MNI_ICBM_2009b_NLIN_ASYM" to "MNI152NLin2009bAsym" ...');
    copyfile(fullfile(LeadRoot, 'templates', 'space', 'MNI_ICBM_2009b_NLIN_ASYM'), fullfile(LeadRoot, 'templates', 'space', 'MNI152NLin2009bAsym'));

    if isfile([ea_space, 'IITmean_tensor.nii.gz'])
        disp('Rename IIT Mean Tensor ...');
        movefile([ea_space, 'IITmean_tensor.nii.gz'], [ea_space, 'IITmeanTensor.nii.gz'])
    end

    if isfile([ea_space, 'IITmean_tensor_Norm_mapping.nii.gz'])
        disp('Rename scaled IIT Mean Tensor ...');
        movefile([ea_space, 'IITmean_tensor_Norm_mapping.nii.gz'], [ea_space, 'IITmeanTensorNormMapping.nii.gz'])
    end

    try
        disp('Discard local changes in LeadDBS repository ...')
        system(['git -C ', LeadRoot, ' checkout .']);
    catch
        warning('off', 'backtrace');
        warning('Failed to discard local changes! Please run "git checkout ." in LeadDBS folder in your Terminal.');
        warning('on', 'backtrace');
    end
end

if isfile([LeadRoot, 'common', filesep, 'ea_recentgroups.mat'])
    disp('Backup recent groups ...');
    movefile([LeadRoot, 'common', filesep, 'ea_recentgroups.mat'], [LeadRoot, 'common', filesep, 'ea_recentgroups.mat.dev'])
end

if isfile([LeadRoot, 'common', filesep, 'ea_recentgroups.mat.bids'])
    disp('Restore recent groups ...');
    movefile([LeadRoot, 'common', filesep, 'ea_recentgroups.mat.bids'], [LeadRoot, 'common', filesep, 'ea_recentgroups.mat'])
end

if isfile([LeadRoot, 'common', filesep, 'ea_recentpatients.mat'])
    disp('Backup recent patients ...');
    movefile([LeadRoot, 'common', filesep, 'ea_recentpatients.mat'], [LeadRoot, 'common', filesep, 'ea_recentpatients.mat.dev'])
end

if isfile([LeadRoot, 'common', filesep, 'ea_recentpatients.mat.bids'])
    disp('Restore recent patients ...');
    movefile([LeadRoot, 'common', filesep, 'ea_recentpatients.mat.bids'], [LeadRoot, 'common', filesep, 'ea_recentpatients.mat'])
end

if isfile([LeadRoot, 'ea_ui.mat'])
    disp('Backup ea_ui.mat ...');
    movefile([LeadRoot, 'ea_ui.mat'], [LeadRoot, 'ea_ui.mat.dev'])
end

if isfile([LeadRoot, 'ea_ui.mat.bids'])
    disp('Restore ea_ui.mat ...');
    movefile([LeadRoot, 'ea_ui.mat.bids'], [LeadRoot, 'ea_ui.mat'])
end

if isfile([ea_gethome, '.ea_prefs.m'])
    disp('Backup .ea_prefs.m ...');
    movefile([ea_gethome, '.ea_prefs.m'], [ea_gethome, '.ea_prefs.m.dev'])
end

if isfile([ea_gethome, '.ea_prefs.m.bids'])
    disp('Restore .ea_prefs.m ...');
    movefile([ea_gethome, '.ea_prefs.m.bids'], [ea_gethome, '.ea_prefs.m'])
end

if isfile([ea_gethome, '.ea_prefs.mat'])
    disp('Backup .ea_prefs.mat ...');
    movefile([ea_gethome, '.ea_prefs.mat'], [ea_gethome, '.ea_prefs.mat.dev'])
end

if isfile([ea_gethome, '.ea_prefs.mat.bids'])
    disp('Restore .ea_prefs.mat ...');
    movefile([ea_gethome, '.ea_prefs.mat.bids'], [ea_gethome, '.ea_prefs.mat'])
end

if isfile([ea_gethome, '.ea_prefs.json'])
    disp('Backup .ea_prefs.json ...');
    movefile([ea_gethome, '.ea_prefs.json'], [ea_gethome, '.ea_prefs.json.dev'])
end

if isfile([ea_gethome, '.ea_prefs.json.bids'])
    disp('Restore .ea_prefs.json ...');
    movefile([ea_gethome, '.ea_prefs.json.bids'], [ea_gethome, '.ea_prefs.json'])
end
