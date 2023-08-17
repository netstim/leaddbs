function ea_leador_create_or_scene(~,~,handles)

if ischar(handles)
    uipatdirs={handles};
else
    uipatdirs=getappdata(handles.leadfigure,'uipatdir');
end

slicer_netstim_modules = {};
slicer_netstim_modules{end+1} = fullfile(ea_getearoot, 'ext_libs', 'SlicerNetstim', 'StereotacticPlan');
slicer_netstim_modules{end+1} = fullfile(ea_getearoot, 'ext_libs', 'SlicerNetstim', 'NetstimPreferences');
slicer_netstim_modules{end+1} = fullfile(ea_getearoot, 'ext_libs', 'SlicerNetstim', 'ImportAtlas');

s4l = ea_slicer_for_lead;
if ~s4l.is_installed_and_up_to_date()
    s4l.install();
end


for i = 1:length(uipatdirs)
    
    pt_options = ea_getptopts(uipatdirs{i});
    subjectAnatFiles = cellfun(@(x) pt_options.subj.coreg.anat.preop.(x), fieldnames(pt_options.subj.coreg.anat.preop), 'uni', 0);
    subjectDICOMFolder = strrep(uipatdirs{i}, [filesep 'derivatives' filesep 'leaddbs' filesep], [filesep 'sourcedata' filesep]);
    MNIToAnchorNativeFile = dir([pt_options.subj.norm.transform.inverseBaseName '*']);
    MNIToAnchorNativeFile = fullfile(MNIToAnchorNativeFile(1).folder, MNIToAnchorNativeFile(1).name);
    atlasPath = fullfile(ea_space, 'atlases', handles.atlassetpopup.String{handles.atlassetpopup.Value}, 'atlas_index.mat');
    
    script_source = [mfilename('fullpath') '.py'];
    script_to_run = fullfile(uipatdirs{i}, 'leador', [mfilename '.py']);
    copyfile(script_source, script_to_run)
    
    lines_to_append = [strcat("subjectAnatFiles=[r'", strjoin(subjectAnatFiles,"r','"), "']"),...
                       strcat("subjectDICOMFolder=r'", subjectDICOMFolder, "'"),...
                       strcat("MNIToAnchorNativeFile=r'", MNIToAnchorNativeFile, "'"),...
                       strcat("atlasPath=r'", atlasPath, "'")];
    
    script_read = fileread(script_to_run);
    script_read = [char(strjoin(lines_to_append,newline)), newline, script_read];
    fileID = fopen(script_to_run, 'w');
    fwrite(fileID, script_read, 'char');
    fclose(fileID);

    command = [' --no-splash'...
           ' --ignore-slicerrc'...
           ' --no-main-window'...
           ' --additional-module-paths "' char(strjoin(string(slicer_netstim_modules),'" "')) '"'...
           ' --python-script "' script_to_run '"'];

    save_log = [' >> "' fullfile(uipatdirs{i}, 'leador', ['sub-' pt_options.subj.subjId '_desc-createORScene' char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss.SSS')) '.txt']) '"'];

    s4l.run([command save_log]);
    delete(script_to_run)
    if ~isfile(fullfile(uipatdirs{i}, 'leador', ['sub-' pt_options.subj.subjId '_desc-ORScene.mrb']))
        warning(['Apparently something went wrong and the Scene file was not created. See the log file in the leador subject folder: ' uipatdirs{i}]);
    else
        display(['Finished sub: ' uipatdirs{i}]);
    end

end
