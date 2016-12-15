function ea_setpath

if ~isdeployed
    addpath(genpath(earoot));
    rmpath(genpath([earoot,'.git']));
    rmpath(genpath([earoot,'release']));
end