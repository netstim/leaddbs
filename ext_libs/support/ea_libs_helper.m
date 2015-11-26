function ea_libs_helper(libpath)

% add shared libraries to search path

if nargin < 1
    libpath = fileparts(mfilename('fullpath'));
end

if ismac
    setenv('DYLD_LIBRARY_PATH',[libpath,':',getenv('DYLD_LIBRARY_PATH')]);
elseif isunix
    setenv('LD_LIBRARY_PATH',[libpath,':',getenv('LD_LIBRARY_PATH')]);
elseif ispc
    setenv('PATH', [libpath, ';', getenv('PATH')]);
end
