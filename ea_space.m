function [path]=ea_space(options,cmd)

spacename='mni_icbm2009b'; % as default for now

switch cmd
    case 'space'
        
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep];
        
    case 'atlases'
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep,'atlases',filesep];
    case {'subcortical','schoenecker'} 
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep,'subcortical',filesep];
    case 'cortex'
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep,'atlases',filesep];
        
    case 'labeling'
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep,'atlases',filesep];

    case 'dartel'
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep,'atlases',filesep];

        
end

if ~exist(path,'file')
    ea_error('This functionality seems not to be compatible with the space you are working in. Please consider using the default ICBM 2009b nonlinear asymmetric space for this procedure.');
end