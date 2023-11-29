function ea_switch2dev

LeadRoot = ea_getearoot;

if ~isfolder(fullfile(LeadRoot, 'templates', 'space', 'MNI152NLin2009bAsym', 'atlases'))
    disp('Copy "MNI_ICBM_2009b_NLIN_ASYM" to "MNI152NLin2009bAsym" ...');
    copyfile(fullfile(LeadRoot, 'templates', 'space', 'MNI_ICBM_2009b_NLIN_ASYM', '*'), fullfile(LeadRoot, 'templates', 'space', 'MNI152NLin2009bAsym'));

    newSpace = [fileparts(fileparts(ea_space)), filesep, 'MNI152NLin2009bAsym', filesep];
    if isfile([newSpace, 'IITmean_tensor.nii.gz'])
        disp('Rename IIT Mean Tensor ...');
        movefile([newSpace, 'IITmean_tensor.nii.gz'], [newSpace, 'IITMeanTensor.nii.gz'])
    end

    if isfile([newSpace, 'IITmean_tensor_Norm_mapping.nii.gz'])
        disp('Rename scaled IIT Mean Tensor ...');
        movefile([newSpace, 'IITmean_tensor_Norm_mapping.nii.gz'], [newSpace, 'IITMeanTensor_NormMapping.nii.gz'])
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

if isfile([ea_prefsdir, filesep, 'ea_recentgroups.mat'])
    disp('Backup recent groups from classic branch ...');
    movefile([ea_prefsdir, filesep, 'ea_recentgroups.mat'], [ea_prefsdir, filesep, 'ea_recentgroups.mat.classic'])
end

if isfile([ea_prefsdir, filesep, 'ea_recentgroups.mat.dev'])
    disp('Restore recent groups from develop branch  ...');
    movefile([ea_prefsdir, filesep, 'ea_recentgroups.mat.dev'], [ea_prefsdir, filesep, 'ea_recentgroups.mat'])
end

if isfile([ea_prefsdir, filesep, 'ea_recentpatients.mat'])
    disp('Backup recent patients from classic branch  ...');
    movefile([ea_prefsdir, filesep, 'ea_recentpatients.mat'], [ea_prefsdir, filesep, 'ea_recentpatients.mat.classic'])
end

if isfile([ea_prefsdir, filesep, 'ea_recentpatients.mat.dev'])
    disp('Restore recent patients from develop branch  ...');
    movefile([ea_prefsdir, filesep, 'ea_recentpatients.mat.dev'], [ea_prefsdir, filesep, 'ea_recentpatients.mat'])
end

if isfile([ea_prefsdir, filesep, 'ea_ui.mat'])
    disp('Backup ea_ui.mat from classic branch  ...');
    movefile([ea_prefsdir, filesep, 'ea_ui.mat'], [ea_prefsdir, filesep, 'ea_ui.mat.classic'])
end

if isfile([ea_prefsdir, filesep, 'ea_ui.mat.dev'])
    disp('Restore ea_ui.mat from develop branch  ...');
    movefile([ea_prefsdir, filesep, 'ea_ui.mat.dev'], [ea_prefsdir, filesep, 'ea_ui.mat'])
end

if isfile(ea_prefspath)
    disp('Backup ea_prefs_user.m from classic branch  ...');
    movefile(ea_prefspath, ea_prefspath('.m.classic'))
end

if isfile(ea_prefspath('.m.dev'))
    disp('Restore ea_prefs_user.m from develop branch  ...');
    movefile(ea_prefspath('.m.dev'), ea_prefspath)
end

if isfile(ea_prefspath('mat'))
    disp('Backup ea_prefs_user.mat from classic branch ...');
    movefile(ea_prefspath('mat'), ea_prefspath('.mat.classic'))
end

if isfile(ea_prefspath('.mat.dev'))
    disp('Restore ea_prefs_user.mat from develop branch  ...');
    movefile(ea_prefspath('.mat.dev'), ea_prefspath('mat'))
    load(ea_prefspath('mat'), 'machine');
    machine.d2.backdrop = 'MNI152NLin2009bAsym T1 (Fonov)';
    machine.togglestates.template = 'MNI152NLin2009bAsym T1 (Fonov)';
    save(ea_prefspath('mat'), 'machine');
end

if isfile(ea_prefspath('json'))
    disp('Backup ea_prefs_user.json from classic branch  ...');
    movefile(ea_prefspath('json'), ea_prefspath('.json.classic'))
end

if isfile(ea_prefspath('.json.dev'))
    disp('Restore ea_prefs_user.json from develop branch  ...');
    movefile(ea_prefspath('.json.dev'), ea_prefspath('json'))
end

disp('Switch LeadDBS branch to develop ...')
system(['git -C ', LeadRoot, ' stash']);
system(['git -C ', LeadRoot, ' checkout develop']);

ea_setpath;
rehash toolboxcache;
disp('LeadDBS search path updated.');
