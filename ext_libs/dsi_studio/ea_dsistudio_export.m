function exportedFile = ea_dsistudio_export(fibgz, data)
% Warpper function to export data (qa/nqa/iso...) from .fib.gz file

basedir = [ea_getearoot, 'ext_libs',filesep,'dsi_studio',filesep];
DSISTUDIO = ea_getExec([basedir, 'dsi_studio'], escapePath = 1);


exportcmd = [DSISTUDIO, ' --action=exp --source=', ea_path_helper(GetFullPath(fibgz)), ' --export=', data];

[~, cmdout] = ea_runcmd(exportcmd);

disp(cmdout);

exportedFile = strip(regexp(cmdout, '(?<=saved to ).*', 'match', 'once'));

if isempty(exportedFile)
    error('"%s" not found in file.', data);
end
