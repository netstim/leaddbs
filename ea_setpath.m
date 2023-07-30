function ea_setpath

% Add lead dir and subdirs
addpath(genpath(ea_getearoot));
rmpath(genpath([ea_getearoot,'.git']));
rmpath(genpath([ea_getearoot,'release']));
rmpath(genpath([ea_getearoot,'ext_libs',filesep,'mambaforge']));
rmpath(genpath([ea_getearoot,'ext_libs',filesep,'fastsurfer',filesep,'upstream']));
rmpath(genpath([ea_getearoot,'ext_libs',filesep,'surfice',filesep,'surfice.app']));
rmpath(genpath([ea_getearoot,'ext_libs',filesep,'SlicerNetstim']));
rmpath(genpath([ea_getearoot,'ext_libs',filesep,'SlicerForLeadDBS']));

% Add SPM dir
% addpath(spm('dir'))

savepath;
