function ea_libs_helper(libpath)

% add shared libraries to search path

if nargin < 1
    libpath = fileparts(mfilename('fullpath'));
end

if ispc
    envname = 'PATH';
elseif ismac
    envname = 'DYLD_LIBRARY_PATH';
elseif isunix
    envname = 'LD_LIBRARY_PATH';
end

env = getenv(envname);

if ~contains(env, libpath)
    setenv(envname, [libpath, ';', env]);
end
