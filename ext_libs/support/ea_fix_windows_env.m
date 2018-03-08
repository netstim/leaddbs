function ea_fix_windows_env
% On Windows, '%MATLABROOT%\bin\win64' might be missing in the
% 'PATH' environment variable. This function will fix it when needed.
% The path is neccesary to run some binaries.

if ispc
    env = getenv('PATH');
    libpath = [matlabroot, filesep, 'bin', filesep, 'win64'];
    if isempty(strfind(env, libpath))
        setenv('PATH', [env, ';', libpath]);
    end
end
