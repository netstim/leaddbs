function [] = ea_runwarpdrive(options)

%
% Get all subjects and only run warpdrive in the last one
%

global WARPDRIVE_SUBS;
if isempty(WARPDRIVE_SUBS)
    WARPDRIVE_SUBS = options.subj;
else
    WARPDRIVE_SUBS(end+1) = options.subj;
end
if options.pat < length(options.uipatdirs)
    return
end
warpdrive_subs = WARPDRIVE_SUBS;
clear global WARPDRIVE_SUBS;

%
% Check which subjects should run
%

remove_pts = false(length(warpdrive_subs),1);
if ~ (isfield(options, 'overwriteapproved') && options.overwriteapproved)
    for i = 1:length(warpdrive_subs)
        if isfile(warpdrive_subs(i).norm.log.method)
            approved_load = loadjson(warpdrive_subs(i).norm.log.method);
            remove_pts(i) = ~isfield(approved_load,'approval') | (approved_load.approval == 0);
            if ~contains(approved_load.method, 'ANTs')
                remove_pts(i) = 1;
                disp([warpdrive_subs(i).subjId ' was normalized using ' approved_load.method '. Use ANTs in order to run warpdrive.']);
            end
        end
    end
end

warpdrive_subs(remove_pts) = [];
if isempty(warpdrive_subs) % all subjects approved
    return
end

%
% Check slicer install
%

slicer_path = ea_runslicer(options, 5);


%
% Check warpdrive install
%

if ~isfield(options.prefs, 'slicer_netstim_path') || strcmp(options.prefs.slicer_netstim_path,'') || ~isfolder(options.prefs.slicer_netstim_path)
    cmd = ['"' slicer_path '"' ...
        ' --no-splash' ...
        ' --no-main-window' ...
        ' --ignore-slicerrc'...
        ' --python-script ' fullfile(ea_getearoot, 'support_scripts', 'install_warpdrive.py')];
    disp('Installing WarpDrive...')
    [status,result] = system(cmd);
    info = regexp(result,'.*os:(?<os>\w+).*rev:(?<rev>\d+).*extensionsInstallPath:(?<extensionsInstallPath>[^\n]+).*slicerMajorVersion:(?<slicerMajorVersion>\d+).*','names');
    if info.slicerMajorVersion < 5
        fprintf('Installed Slicer version is %d. This might couse some issues. We recommend to update to Slicer 5.\n <a href = "%s">%s</a>\n', info.slicerMajorVersion, 'https://download.slicer.org/', 'https://download.slicer.org/');
    end
    switch status
        case 0
            d = dir(fullfile(info.extensionsInstallPath,'**','SlicerNetstim'));
            slicer_netstim_path = d(1).folder;
            ea_injectprefstring('slicer_netstim_path', slicer_netstim_path);
            options.prefs.slicer_netstim_path = slicer_netstim_path;
            disp('WarpDrive installed!');
        case 1
            warndlg('WarpDrive could not be installed automatically, probably because of a network related issue. See the command window output for options on how to install manually.', 'Warning');
            if regexp(result,'.*Execution of PAC script at.*proxy.*') && startsWith(computer('arch'),'mac')
                slicer_opt_page = 'https://github.com/simonoxen/SlicerForMacWithProxy/releases/tag/v5.0.2';
                fprintf('Looks like you are using MacOS under a proxy. You can download a Slicer version that fixes a proxy-related bug here:\n<a href = "%s">%s</a>\n',slicer_opt_page,slicer_opt_page);
            else
                download_page = ['https://extensions.slicer.org/view/SlicerNetstim/' info.rev '/' info.os];
                extenstions_manager_page = 'https://slicer.readthedocs.io/en/latest/user_guide/extensions_manager.html#install-downloaded-extension-packages';
                fprintf(['You can manually download the extension from here:\n<a href = "%s">%s</a>\n' ...
                    'And intall it through the extensions manager:\n<a href = "%s">%s</a>\n'], download_page, download_page, extenstions_manager_page, extenstions_manager_page);
            end
            disp('If this error persists contact us via our user slack channel');
            return
    end
end

%
% Get the subject's necesary information to send to warpdrive
%

info_struct = struct;

for i = 1:length(warpdrive_subs)
    
    info_struct(i).id = warpdrive_subs(i).subjId;
    info_struct(i).warpdrive_path = warpdrive_subs(i).warpdriveDir;
    info_struct(i).normlog_file = warpdrive_subs(i).norm.log.method;
    info_struct(i).anat_files = warpdrive_subs(i).coreg.anat.preop;
    
    d = dir([warpdrive_subs(i).norm.transform.forwardBaseName, 'ants*']);
    if endsWith(d(1).name, '.h5')
        update_ants_transforms(warpdrive_subs(i));
    end

    info_struct(i).forward_transform = [warpdrive_subs(i).norm.transform.forwardBaseName, 'ants.nii.gz'];
    info_struct(i).inverse_transform = [warpdrive_subs(i).norm.transform.inverseBaseName, 'ants.nii.gz'];
   
end

tmp_file = fullfile(options.subj.logDir, '.warpdrive_tmp.json');
fid = fopen(tmp_file, 'w');
fprintf(fid,'%s', jsonencode(info_struct));
fclose(fid);

%
% Set up commands and run warpdrive
%

slicer_netstim_modules = dir(fullfile(options.prefs.slicer_netstim_path,'**','cli-modules'));
slicer_netstim_modules = unique({slicer_netstim_modules.folder}');
slicer_netstim_modules{end+1} = fullfile(ea_getearoot, 'ext_libs', 'SlicerNetstim', 'ImportAtlas');
slicer_netstim_modules{end+1} = fullfile(ea_getearoot, 'ext_libs', 'SlicerNetstim', 'NetstimPreferences');
slicer_netstim_modules{end+1} = fullfile(ea_getearoot, 'ext_libs', 'SlicerNetstim', 'WarpDrive');

python_commands = [strcat("slicer.app.settings().setValue('NetstimPreferences/leadDBSPath',r'", remove_last_filesep(ea_getearoot), "')");...
                    "slicer.app.settings().setValue('MainWindow/DontShowDisclaimerMessage','1024')";...
                    "slicer.app.settings().setValue('MainWindow/DontConfirmExit','1024')";...
                    "slicer.app.settings().setValue('MainWindow/DontConfirmRestart','1024')";...
                    "import json, os, sys";...
                    "list(map(lambda p: sys.path.insert(0,sys.path.pop(sys.path.index(p))), [p for p in sys.path if p.find('leaddbs')>0]))";...
                    "import WarpDrive";...
                    strcat("WarpDrive.WarpDriveLogic().getParameterNode().SetParameter('LeadSubjects',json.dumps(json.load(open(r'", tmp_file, "'))))");...
                    strcat("WarpDrive.WarpDriveLogic().getParameterNode().SetParameter('MNIPath',r'", remove_last_filesep(ea_space), "')");...
                    strcat("os.remove(r'", tmp_file, "')");...
                    "slicer.util.selectModule('WarpDrive')"];
                   
command = ['"' slicer_path '"' ...
           ' --no-splash'...
           ' --disable-settings'...
           ' --ignore-slicerrc'...
           ' --additional-module-paths "' char(strjoin(string(slicer_netstim_modules),'" "')) '"'...
           ' --python-code "' char(strjoin(python_commands,";")) '"'];

save_log = [' >> "' fullfile(warpdrive_subs(1).logDir, ['warpdrive_' char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss.SSS')) '.txt']) '"'];
       
system([command save_log ' &']); % with & return control to Matlab
disp('Running WarpDrive in Slicer');

end

function out = remove_last_filesep(filepath)
    if strcmp(filepath(end), filesep)
        filepath = filepath(1:end-1);
    end
    out = filepath;
end


function [] = update_ants_transforms(subj)
    fprintf('Updating transform from .h5 to .nii.gz for subject: %s\n', subj.subjId);

    if ispc
        ext = 'exe';
    else
        ext = computer('arch');
    end
    ants_apply = fullfile(ea_getearoot, 'ext_libs', 'ANTs', ['antsApplyTransforms.' ext]);
    
    transforms_base = {subj.norm.transform.forwardBaseName, subj.norm.transform.inverseBaseName};
    references = {fullfile(ea_space, 't1.nii'), subj.coreg.anat.preop.(subj.AnchorModality)};
    
    for j = 1:length(transforms_base)
        
        t = [transforms_base{j} 'ants.h5'];
        o = ['[' transforms_base{j} 'ants.nii.gz,1]'];
        r = references{j};
        
        cmd = [ants_apply, ' -r ' r ' -o ' o ' -t ' t ' -v 1 --float'];
        if ~ispc
            system(['bash -c "', cmd, '"']);
        else
            system(cmd);
        end
        
        delete(t);
        
    end
end