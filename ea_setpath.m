function ea_setpath

if ~isdeployed
    addpath(genpath(ea_getearoot));
    rmpath(genpath([ea_getearoot,'.git']));
    rmpath(genpath([ea_getearoot,'release']));
end