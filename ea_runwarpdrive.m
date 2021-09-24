function [] = ea_runwarpdrive(options)

if options.pat > 1 % run all subjects at once with the first one
    return
end

do_pts = true(length(options.uipatdirs),1);

if ~ (isfield(options, 'overwriteapproved') && options.overwriteapproved)
    
    for i = 1:length(options.uipatdirs)
        approved_file = dir(fullfile(options.uipatdirs{i},'normalization','log','*desc-normmethod.json'));
        if ~isempty(approved_file)
            approved_file = fullfile(approved_file(1).folder, approved_file(1).name);
            approved_load = loadjson(approved_file);
            do_pts = ~isfield(approved_load,'approval') | (approved_load.approval == 0);
        end
    end
    
end


if ~any(do_pts) % all subjects approved
    return
end

do_pts_dirs = options.uipatdirs(do_pts);

% get slicer / slicer custom executable
d = dir(fullfile(ea_getearoot,'ext_libs','SlicerCustom*'));
if isempty(d)
    slicer_path = ea_runslicer(options, 5);
    save_log = '';
else
    d = d([d.isdir]); % keep directories only
    if ismac
        slicer_path = fullfile(d(1).folder,d(1).name,'SlicerCustom.app','Contents','MacOS','SlicerCustom');
    elseif isunix
        % TODO
    elseif ispc
        % TODO
    end
    % set up to save log file
    log_path = fullfile(d(1).folder, d(1).name, 'log');
    if ~isfolder(log_path)
        mkdir(log_path)
    end
    save_log = [' >> "' fullfile(log_path, [char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss.SSS')) '.txt']) '"'];
end

% aditional modules
addition_module_paths = {fullfile(ea_getearoot,'ext_libs','SlicerNetstim','ImportAtlas'),...
                         fullfile(ea_getearoot,'ext_libs','SlicerNetstim','WarpDrive')}; 

command = ['"' slicer_path '"' ...
           ' --no-splash'...
           ' --ignore-slicerrc'...
           ' --additional-module-paths "' strjoin(addition_module_paths,'" "') '"' ...        % SlicerNetstim modules 
           ' --python-code "slicer.util.selectModule(''WarpDrive'')" ' ...  % Change to WarpDrive module
           ' "' strjoin([ea_getearoot; do_pts_dirs],'" "') '"'];               % Additional args with leadroot and pts dir
       
system([command save_log ' &']); % with & return control to Matlab
disp('Running WarpDrive in Slicer');

end