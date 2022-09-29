function ea_setpath

% Add lead dir and subdirs
addpath(genpath(ea_getearoot));
rmpath(genpath([ea_getearoot,'.git']));
rmpath(genpath([ea_getearoot,'release']));
rmpath(genpath([ea_getearoot,'ext_libs',filesep,'mambaforge']));

% Add SPM dir
% addpath(spm('dir'))

