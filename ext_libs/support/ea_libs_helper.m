function ea_libs_helper(libpath, setpath)

% add/remove shared libraries to search path

% Use current folder as libpath by default
if nargin < 1 || isempty(libpath)
    %load prefs for platform specific handling
    prefs = ea_prefs;

    % Check and load runtime libs when needed
    arch = ea_getarch;
    if eval(['prefs.platform.', arch, '.load_shipped_runtime'])
        libpath = fullfile(fileparts(mfilename('fullpath')), arch);
    else
        % Nothing to do, exit
        return;
    end
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
            if ispc()
                setenv(envname, [libpath, ';', env]); % For pc is ";"
            else
                setenv(envname, [libpath, ':', env]); % For unix (including mac, is ":")
            end
        end
    case {0, 'unset'}
        env = regexprep(env, [libpath,'[;:]'], '');
        setenv(envname, env);
end
