function [] = ea_runwarpdrive(options)

if options.pat > 1 % run all subjects at once with the first one
    return
end

do_pts = true(length(options.uipatdirs),1);

if ~ (isfield(options, 'overwriteapproved') && options.overwriteapproved)
    
    for i = 1:length(options.uipatdirs)
        approved_file = fullfile(options.uipatdirs{i},'ea_coreg_approved.mat');
        if isfile(approved_file)
            approved_load = load(approved_file);
            if isfield(approved_load,'glanat') && approved_load.glanat
                do_pts(i) = false;
            end
        end
    end
    
end


if ~any(do_pts) % all subjects approved
    return
end

do_pts_dirs = options.uipatdirs(do_pts);

slicer_path = ea_runslicer(options,5);
d = dir(fullfile(ea_getearoot,'ext_libs','SlicerNetstim','*','CMakeLists.txt')); % aditional modules

command = [slicer_path ' --no-splash --additional-module-paths ' strjoin({d.folder},' ') ' --python-code "slicer.util.selectModule(''WarpDrive'')" ' ea_getearoot ' "' strjoin(do_pts_dirs,'" "') '"'];
system([command ' &']); % with & return control to Matlab
disp('Running WarpDrive in Slicer');

end