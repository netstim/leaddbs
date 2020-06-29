function mask=ea_getmask(varargin)
% shortcut to quickly get the commonly used brainmask used in 2x2x2 within
% predict suite.
if nargin
    switch varargin{1}
        case {'gray','grey','','standard','grey_2'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.mask;
        case {'gray_hd','grey_hd','hd','standard_hd','grey_1','gray_1'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.mask_hd;
        case {'brain','brain_2'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.brainmask;
        case {'brain_hd','brain_1'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.brainmask_hd;
        case {'cortex','cortex_2'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.cortexmask;
        case {'cortexcb','cortexcb_2'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.cortexcbmask;
        case {'cortex_hd','cortex_1'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.cortexmask_hd;
        case {'cortexcb_hd','cortexcb_1'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.cortexcbmask_hd;
        case {'cortex_5'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.cortexmask_hd5;
        case {'brain_5'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.brainmask_hd5;
        case {'grey_5','gray_5'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.mask_hd5;
        case {'cb','cb_2'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.cbmask;
        case {'cb_hd','cb_1'}
            load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
            mask=modeldata.cbmask_hd;
    end
else
    load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']);
    mask=modeldata.mask;
end