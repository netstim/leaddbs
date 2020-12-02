function ea_libs_helper(libpath, setpath)

% add/remove shared libraries to search path

% Use current folder as libpath by default
if nargin < 1 || isempty(libpath)
    libpath = fileparts(mfilename('fullpath'));
end

% set path by default
if nargin < 2
    setpath = 1;
end

if ispc
    envname = 'PATH';
elseif ismac
    envname = 'DYLD_LIBRARY_PATH';
elseif isunix
    envname = 'LD_LIBRARY_PATH';
end

env = getenv(envname);

switch setpath
    case {1, 'set'}
        if ~contains(env, libpath)
            setenv(envname, [libpath, ';', env]);
        end
    case {0, 'unset'}
        env = regexprep(env, [libpath,'[;:]'], '');
        setenv(envname, env);
end
