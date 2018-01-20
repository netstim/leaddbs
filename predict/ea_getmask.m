function mask=ea_getmask
% shortcut to quickly get the commonly used brainmask used in 2x2x2 within
% predict suite.
load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);

mask=modeldata.mask;