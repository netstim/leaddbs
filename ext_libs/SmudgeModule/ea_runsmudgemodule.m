function [] = ea_runsmudgemodule(options)


slicer_path = ea_runslicer(options,5);
module_path = fullfile(ea_getearoot,'ext_libs','SmudgeModule');
module_script = fullfile(module_path,'SmudgeModule.py');
command = [slicer_path ' --no-main-window --additional-module-paths ' module_path ' --python-script ' module_script ' ' ea_getearoot ' ' strjoin(options.uipatdirs,' ')];
system(command);