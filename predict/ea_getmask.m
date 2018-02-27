function mask=ea_getmask(varargin)
% shortcut to quickly get the commonly used brainmask used in 2x2x2 within
% predict suite.
if nargin
    switch varargin{1}
        case {'grey','','standard'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.mask;
        case 'brain'
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.brainmask;
        case 'cortex'
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.cortexmask;
    end
else
    load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
    mask=modeldata.mask;
end