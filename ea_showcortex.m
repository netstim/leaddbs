function [Hp,cortex] = ea_showcortex(varargin)
% This function shows a cortical reconstruction in the 3D-Scene viewer. It
% reads in cortex.mat found in the eAuto_root/templates/cortex folder, or 
% in the patientdirectory on button press (Cortical Reconstruction
% Visualization) in the Lead-Anatomy scene. To view patient specific 
% cortical reconstructions, load them to the patient directory using
% ea_importfs via the "Import FS" button in the main Lead DBS GUI.
% __________________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh (UPMC), Brain Modulation Lab
% Ari Kappel

resultfig=varargin{1};
if nargin==2
    options=varargin{2};
end
if nargin>2
    elstruct=varargin{2};
    options=varargin{3};
end

% Initial Opening
if ~isfield(getappdata(resultfig),'cortex')
    color = options.prefs.d3.cortexcolor; % default color is gray
    alpha = options.prefs.d3.cortexalpha; % default alph is 0.333
else
    color = options.prefs.d3.cortexcolor;
    appdata = getappdata(resultfig);
    appdata = getappdata(appdata.awin);
    alpha = str2double(appdata.UsedByGUIData_m.cortexalpha.String);
    clear appdata
end

    
nm=[0:2]; % native and mni
try
    nmind=[options.atl.pt,options.atl.can,options.atl.ptnative]; % which shall be performed?
catch
    nmind=[0 1 0];
end
nm=nm(logical(nmind)); % select which shall be performed.

% switch between patient and template cortex
slicecontroldata = getappdata(gcf);
if ~isempty(strfind(slicecontroldata.templateused,'Patient'))
    nm = 0;
end

for nativemni=nm % switch between native and mni space.
    
    switch nativemni
        case 0 % patient cortex in mni space
            root=[options.root,options.patientname,filesep];
            adir=[root,''];
            reslice='yes';
        case 1 % template cortex in mni space
            root=[options.earoot,mcr];
            adir=[root,'templates/cortex/'];
            reslice='no';
        case 2 % patient cortex in native space
            root=[options.root,options.patientname,filesep];
            adir=[root,''];
            reslice='no';
    end
end

try
    load([adir,'cortex.mat'])
catch
    disp('Running Import FS...')
    ea_importfs(options)
    load([adir,'cortex.mat'])
end

% % reslice patient cortex to mni space in future release
% % for now always use template cortex in mni space
% if strcmp(reslice,'yes')
%     
% end

% Show cortex
hold on
Hp = patch('vertices',cortex.vert,'faces',cortex.tri(:,[1 3 2]),...
    'facecolor',color,'edgecolor','none','FaceAlpha',alpha,...
    'facelighting', 'gouraud', 'specularstrength', .25);
% camlight('headlight','infinite'); axis equal;
