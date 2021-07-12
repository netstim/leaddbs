function ea_surfice_script(script,hold)
% calls surfice

if ~exist('hold','var')
    hold = 1; % stop other matlab executions - should use if script contains a QUIT element
end

% check OpenGL version
openglInfo = opengl('data');
openglVer = strrep(regexp(openglInfo.Version, '(^[\d.]+)(?=.*)', 'match', 'once'), '.', '');
if numel(openglVer) == 2
    openglVer = [openglVer, '0'];
end
openglVer = str2double(openglVer);

if openglVer >= 330 || ismac
    surfice_exe = 'surfice';
elseif openglVer >= 210
    surfice_exe = 'surficeOld';
else
    ea_error(sprintf('Surf Ice failed to load proper OpenGL!\nFound version: %s, Vendor: %s.', ...
                      openglInfo.Version, openglInfo.Vendor),'Error', dbstack);
end

basedir=[ea_getearoot,'ext_libs',filesep,'surfice',filesep];
if ismac
    surfice = [basedir, surfice_exe, '.app',filesep,'Contents',filesep,'MacOS',filesep,'surfice'];
elseif isunix
    ea_libs_helper([basedir, 'linux']);
    surfice = [basedir, surfice_exe];
elseif ispc
    surfice = ea_path_helper([basedir, surfice_exe, '.exe']);
end

cmd = [surfice, ' -S "', script,'"'];
if ~hold
    cmd = [cmd,' &'];
end
system(cmd);
