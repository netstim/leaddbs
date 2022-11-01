function ea_surfice(script, hold)
% Wrapper for Surf-Ice script call

if ~exist('hold', 'var')
    hold = 1; % Block executions until script QUIT
end

% check OpenGL version
% openglInfo = opengl('data');
% openglVer = strrep(regexp(openglInfo.Version, '(^[\d.]+)(?=.*)', 'match', 'once'), '.', '');
% if numel(openglVer) == 2
%     openglVer = [openglVer, '0'];
% end
% openglVer = str2double(openglVer);

% if openglVer >= 330 || ismac
    SURFICE = 'surfice';
% elseif openglVer >= 210
%     SURFICE = 'surficeOld';
% else
%     ea_error(sprintf('Surf Ice failed to load proper OpenGL!\nFound version: %s, Vendor: %s.', ...
%         openglInfo.Version, openglInfo.Vendor), 'Error', dbstack);
% end

basedir = [fileparts(mfilename('fullpath')), filesep];

% Get binary path
if ismac
    surfice = ea_path_helper([basedir, SURFICE, '.app',filesep,'Contents',filesep,'MacOS',filesep,'surfice']);
elseif isunix
    surfice = ea_path_helper([basedir, SURFICE]);
elseif ispc
    surfice = ea_path_helper([basedir, SURFICE, '.exe']);
end

% Construct script call
cmd = [surfice, ' -S "', script, '"'];

% Redirect stderr to null for Unix
if isunix
    cmd = [cmd, ' 2> /dev/null'];
end

% when hold is false, run the command in background
if ~hold
    if ispc
        cmd = ['start /B ', cmd];
    else
        cmd = [cmd,' &'];
    end
end

% Handle library conflicts for Linux
if isunix && ~ismac
    libPathBackup = getenv('LD_LIBRARY_PATH');
    setenv('LD_LIBRARY_PATH', strrep(libPathBackup, [matlabroot,'/cefclient/sys/os/glnxa64:'],''));
end

% Run Surf-Ice with script
system(cmd);

% Restore original library path
if isunix && ~ismac
    setenv('LD_LIBRARY_PATH', libPathBackup);
end
