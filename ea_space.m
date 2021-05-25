function [path]=ea_space(options,cmd,native)

spacename=ea_getspace; % as default for now
if ~exist('cmd','var')
    cmd='space';
end

if ~exist('native','var')
    native=0; % additional native variable will be used if not working in native space but still wanting to e.g. list atlases from native space. Is not important.
end

if ~exist('options','var')
   options=struct;
end

if ~isfield(options,'native')
    options.native=0;
end

if exist('native','var') % overwrite options.native with native
   options.native=native;
end

switch cmd
    case 'space'
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep];
    case 'atlases'
        if options.native || native
            path=[options.root,options.patientname,filesep,'atlases',filesep];
        else
            path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep,'atlases',filesep];
        end
    case {'subcortical','schoenecker'}
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep,'subcortical',filesep];
    case 'cortex'
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep,'cortex',filesep];

    case 'labeling'
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep,'labeling',filesep];

    case 'dartel'
        path=[ea_getearoot,'templates',filesep,'space',filesep,spacename,filesep,'dartel',filesep];
    case 'connectomes'
        p=ea_prefs;
        path=p.lc.datadir;
end

if ~exist(path,'dir')
    warning('off', 'backtrace');
    warning('Not all atlas data is present. Please click on ''Install'' -> ''Redownload data files'' to download the required files');
    warning('on', 'backtrace');
    %warning('This functionality seems not to be compatible with the space you are working in. Please consider using the default ICBM 2009b nonlinear asymmetric space for this procedure.');
end
