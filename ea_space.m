function [path]=ea_space(options,cmd)

spacename='mni_icbm2009b'; % as default for now

switch cmd
    case 'space'
        
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep];
        
    case 'atlases'
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep,'atlases',filesep];
        
        
    case 'cortex'
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep,'atlases',filesep];
        
    case 'labeling'
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep,'atlases',filesep];

    case 'dartel'
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep,'atlases',filesep];

        
end