function ea_surfice(script,hold)
% calls surfice

if ~exist('hold','var')
    hold=1; % stop other matlab executions - should use if script contains a QUIT element
end
basedir=[ea_getearoot,'ext_libs',filesep,'surfice',filesep];
if ismac
    surfice = [basedir,'mac',filesep, 'surfice.app',filesep,'Contents',filesep,'MacOS',filesep,'surfice'];
elseif isunix
    ea_libs_helper([basedir, 'linux']);
    surfice = [basedir, 'linux',filesep,'surfice'];
elseif ispc
    surfice = [basedir, 'win',filesep,'surfice.exe'];
end
cmd=[surfice,' -S "',script,'"'];
if ~hold
    cmd=[cmd,' &'];
end
system(cmd);
