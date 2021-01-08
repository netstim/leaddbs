function ea_fix_runtimepath
% In some installation, MATLAB run-time library path might be missing in
% the environment variable. This function will fix it when needed.
% The path is neccesary to run some binaries.

runtimePath = {fullfile(matlabroot,'bin',computer('arch'))};

if ispc
    envname = 'PATH';
else
    if ismac
        envname = 'DYLD_LIBRARY_PATH';
    elseif isunix
        envname = 'LD_LIBRARY_PATH';
    end
    runtimePath = [runtimePath; fullfile(matlabroot, 'sys', 'os', computer('arch'))];
end

env = getenv(envname);

for i=1:length(runtimePath)
    if ~contains(env, runtimePath{i})
        if ispc()
            env = [env, ';', runtimePath{i}]; % For pc is ";"
        else
            env = [env, ':', runtimePath{i}]; % For unix (including mac, is ":")
        end
    end
end

setenv(envname, env);
