function mask=ea_getmask(varargin)
% shortcut to quickly get the commonly used brainmask used in 2x2x2 within
% predict suite.
if nargin
    switch varargin{1}
        case {'grey','','standard'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.mask;
        case {'grey_hd','hd','standard_hd'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.mask_hd;
        case 'brain'
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.brainmask;
        case 'brain_hd'
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.brainmask_hd;
        case 'cortex'
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.cortexmask;
        case 'cortex_hd'
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.cortexmask_hd;
    end
else
    load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
    mask=modeldata.mask;
end