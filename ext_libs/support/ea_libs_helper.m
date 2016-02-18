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

if isempty(strfind(getenv(envname),libpath))
    setenv(envname, [libpath, ';', getenv(envname)]);
end
