function ea_libs_helper(libpath, setpath)

% add/remove shared libraries to search path

% Use current folder as libpath by default
if nargin < 1 || isempty(libpath)
    %load prefs for platform specific handling
    prefs = ea_prefs;
    
    libpath = fileparts(mfilename('fullpath'));
    %load system specific needded libs
    if ispc
        if prefs.platform.pcwin.load_shipped_libstdcpp6
            libpath = fullfile(libpath,'pcwin');
        else
            %nothing to do, exit
            return;
        end
    elseif ismac
        if prefs.platform.mac.load_shipped_libstdcpp6
            libpath = fullfile(libpath,'mac');
        else
            %nothing to do, exit
            return;
        end
    elseif isunix
        if prefs.platform.unix.load_shipped_libstdcpp6
            libpath = fullfile(libpath,'unix');
        else
            %nothing to do, exit
            return;
        end
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
                setenv(envname, [libpath, ';', env]);%for pc is ";"
            else
                setenv(envname, [libpath, ':', env]);%for unix (including mac, is ":")
            end
        end
    case {0, 'unset'}
        env = regexprep(env, [libpath,'[;:]'], '');
        setenv(envname, env);
end
