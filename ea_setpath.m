function ea_setpath

% Add lead dir and subdirs
addpath(genpath(ea_getearoot));
rmpath(genpath([ea_getearoot,'.git']));
rmpath(genpath([ea_getearoot,'release']));

% Add SPM dir
% addpath(spm('dir'))

