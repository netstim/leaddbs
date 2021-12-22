function exportedFile = ea_dsistudio_export(fibgz, data)
% Warpper function to export data (qa/nqa/iso...) from .fib.gz file

basedir = [ea_getearoot, 'ext_libs',filesep,'dsi_studio',filesep];
if ispc
    DSISTUDIO = ea_path_helper([basedir, 'dsi_studio.exe']);
else
    DSISTUDIO = [basedir, 'dsi_studio.', computer('arch')];
end

exportcmd = [DSISTUDIO, ' --action=exp --source=', ea_path_helper(GetFullPath(fibgz)), ' --export=', data];

if ~ispc
    [~, cmdout] = system(['bash -c "', exportcmd, '"']);
else
    [~, cmdout] = system(exportcmd);
end

disp(cmdout);

exportedFile = strip(regexp(cmdout, '(?<=saved to ).*', 'match', 'once'));

if isempty(exportedFile)
    error('"%s" not found in file.', data);
end
