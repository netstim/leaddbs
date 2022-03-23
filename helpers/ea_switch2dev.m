function ea_switch2dev

LeadRoot = ea_getearoot;

if ~isfolder(fullfile(LeadRoot, 'templates', 'space', 'MNI152NLin2009bAsym', 'atlases'))
    disp('Copy "MNI_ICBM_2009b_NLIN_ASYM" to "MNI152NLin2009bAsym" ...');
    copyfile(fullfile(LeadRoot, 'templates', 'space', 'MNI_ICBM_2009b_NLIN_ASYM', '*'), fullfile(LeadRoot, 'templates', 'space', 'MNI152NLin2009bAsym'));

    newSpace = [fileparts(fileparts(ea_space)), filesep, 'MNI152NLin2009bAsym', filesep];
    if isfile([newSpace, 'IITmean_tensor.nii.gz'])
        disp('Rename IIT Mean Tensor ...');
        movefile([newSpace, 'IITmean_tensor.nii.gz'], [newSpace, 'IITmeanTensor.nii.gz'])
    end

    if isfile([newSpace, 'IITmean_tensor_Norm_mapping.nii.gz'])
        disp('Rename scaled IIT Mean Tensor ...');
        movefile([newSpace, 'IITmean_tensor_Norm_mapping.nii.gz'], [newSpace, 'IITmeanTensorNormMapping.nii.gz'])
    end

    try
        disp('Discard local changes in LeadDBS repository ...')
        system(['git -C ', LeadRoot, ' checkout .']);
        [~, cmdout] = system(['git -C ', LeadRoot, ' clean -fn']);
        files = regexp(cmdout, '(?<=Would remove )(templates/space/MNI152NLin2009bAsym/[\w/.]+)', 'match')';
        ea_delete(strcat(LeadRoot, files));
    catch
        warning('off', 'backtrace');
        warning('Failed to discard local changes! Please run "git checkout ." in LeadDBS folder in your Terminal.');
        warning('on', 'backtrace');
    end
end

if isfile([LeadRoot, 'common', filesep, 'ea_recentgroups.mat'])
    disp('Backup recent groups from classic branch ...');
    movefile([LeadRoot, 'common', filesep, 'ea_recentgroups.mat'], [LeadRoot, 'common', filesep, 'ea_recentgroups.mat.classic'])
end

if isfile([LeadRoot, 'common', filesep, 'ea_recentgroups.mat.dev'])
    disp('Restore recent groups from develop branch  ...');
    movefile([LeadRoot, 'common', filesep, 'ea_recentgroups.mat.dev'], [LeadRoot, 'common', filesep, 'ea_recentgroups.mat'])
end

if isfile([LeadRoot, 'common', filesep, 'ea_recentpatients.mat'])
    disp('Backup recent patients from classic branch  ...');
    movefile([LeadRoot, 'common', filesep, 'ea_recentpatients.mat'], [LeadRoot, 'common', filesep, 'ea_recentpatients.mat.classic'])
end

if isfile([LeadRoot, 'common', filesep, 'ea_recentpatients.mat.dev'])
    disp('Restore recent patients from develop branch  ...');
    movefile([LeadRoot, 'common', filesep, 'ea_recentpatients.mat.dev'], [LeadRoot, 'common', filesep, 'ea_recentpatients.mat'])
end

if isfile([LeadRoot, 'ea_ui.mat'])
    disp('Backup ea_ui.mat from classic branch  ...');
    movefile([LeadRoot, 'ea_ui.mat'], [LeadRoot, 'ea_ui.mat.classic'])
end

if isfile([LeadRoot, 'ea_ui.mat.dev'])
    disp('Restore ea_ui.mat from develop branch  ...');
    movefile([LeadRoot, 'ea_ui.mat.dev'], [LeadRoot, 'ea_ui.mat'])
end

if isfile([ea_gethome, '.ea_prefs.m'])
    disp('Backup .ea_prefs.m from classic branch  ...');
    movefile([ea_gethome, '.ea_prefs.m'], [ea_gethome, '.ea_prefs.m.classic'])
end

if isfile([ea_gethome, '.ea_prefs.m.dev'])
    disp('Restore .ea_prefs.m from develop branch  ...');
    movefile([ea_gethome, '.ea_prefs.m.dev'], [ea_gethome, '.ea_prefs.m'])
end

if isfile([ea_gethome, '.ea_prefs.mat'])
    disp('Backup .ea_prefs.mat from classic branch ...');
    movefile([ea_gethome, '.ea_prefs.mat'], [ea_gethome, '.ea_prefs.mat.classic'])
end

if isfile([ea_gethome, '.ea_prefs.mat.dev'])
    disp('Restore .ea_prefs.mat from develop branch  ...');
    movefile([ea_gethome, '.ea_prefs.mat.dev'], [ea_gethome, '.ea_prefs.mat'])
end

if isfile([ea_gethome, '.ea_prefs.json'])
    disp('Backup .ea_prefs.json from classic branch  ...');
    movefile([ea_gethome, '.ea_prefs.json'], [ea_gethome, '.ea_prefs.json.classic'])
end

if isfile([ea_gethome, '.ea_prefs.json.dev'])
    disp('Restore .ea_prefs.json from develop branch  ...');
    movefile([ea_gethome, '.ea_prefs.json.dev'], [ea_gethome, '.ea_prefs.json'])
end

disp('Switch LeadDBS branch to develop ...')
system(['git -C ', LeadRoot, ' checkout .']);
system(['git -C ', LeadRoot, ' checkout develop']);

lead path;
rehash toolboxcache;
disp('LeadDBS search path updated.');
