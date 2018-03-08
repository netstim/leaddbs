function ea_libs_helper(libpath)

% add shared libraries to search path

if nargin < 1
    libpath = fileparts(mfilename('fullpath'));
end

if ismac
    envname = 'DYLD_LIBRARY_PATH';
elseif isunix
    envname = 'LD_LIBRARY_PATH';
elseif ispc
    envname = 'PATH';
end

env = getenv(envname);

if isempty(strfind(env, libpath))
    setenv(envname, [libpath, ';', env]);
end
