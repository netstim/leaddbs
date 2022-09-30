function electrode=ea_elspec_smartflow_ngs_nc_06(varargin)
% Creates the FEM-compatible electrode model for SmartFlow Cannula.
% It's based on the mesh generated using SketchUp and tetgen.
% _________________________________________________________________________
% Copyright (C) 2022 Charite University Medicine Berlin
% Ningfei Li

% Set folder
elemodelPath = fileparts(mfilename('fullpath'));
modelFolder = 'SmartFlow_NGS-NC-06';

% Get specification
options.elmodel='SmartFlow Cannula NGS-NC-06';
options = ea_resolve_elspec(options);
elspec = options.elspec;

% Get insulation and contact numbers
numIns = numel(ea_regexpdir([elemodelPath, filesep, modelFolder,filesep, 'Insulations'], '.*\.smesh$'));
numCon = numel(ea_regexpdir([elemodelPath, filesep, modelFolder,filesep, 'Contacts'], '.*\.smesh$'));

%% Import insulations and contacts meshes
for k = 1:numIns
    filename = [elemodelPath, filesep, modelFolder, filesep, 'Insulations', filesep, 'ins', num2str(k), '.1'];
    [node,~,face]=readtetgen(filename);
    electrode.insulation(k).vertices = node;
    electrode.insulation(k).faces = face(:,1:3);
    clear face node filename
end

for k = 1:numCon
    filename = [elemodelPath, filesep, modelFolder, filesep, 'Contacts', filesep, 'con', num2str(k), '.1'];
    [node,~,face]=readtetgen(filename);
    electrode.contacts(k).vertices = node;
    electrode.contacts(k).faces = face(:,1:3);
    clear face node filename
end

%% Contact coordinates and other specifications
electrode.coords_mm(1,:)=[0 0 elspec.tip_length/2];
electrode.coords_mm(2,:)=[0 0 elspec.tip_length+elspec.contact_length/2];

electrode.head_position = [0 0 elspec.tip_length/2]; % dummy value
electrode.tail_position = [0 0 elspec.tip_length+elspec.contact_length/2]; % dummy value
electrode.x_position = [elspec.tip_diameter/2, 0, elspec.tip_length/2]; % dummy value
electrode.y_position = [0, elspec.tip_diameter/2, elspec.tip_length/2]; % dummy value

electrode.electrode_model = options.elmodel;
electrode.isdirected = elspec.isdirected;
electrode.numel = elspec.numel;
electrode.contact_color = elspec.contact_color;
electrode.lead_color = elspec.lead_color;

%% saving electrode struct
save([elemodelPath, filesep, elspec.matfname, '.mat'],'electrode');

%% create and save _vol file
filename = [elemodelPath, filesep, modelFolder, filesep, 'final.1'];
[node,~,face] = readtetgen(filename);
save([elemodelPath, filesep, elspec.matfname, '_vol.mat'],'face','node')
clear node face

%% visualize
if nargin
    vizz=0;
else
    vizz=1;
end

if vizz
    X = eye(4);
    aData = 1;

    figure;
    for ins=1:length(electrode.insulation)
        vs=X*[electrode.insulation(ins).vertices,ones(size(electrode.insulation(ins).vertices,1),1)]';
        electrode.insulation(ins).vertices=vs(1:3,:)';
        elrender=patch('Faces',electrode.insulation(ins).faces,'Vertices',electrode.insulation(ins).vertices);
        ea_specsurf(elrender,elspec.lead_color,aData);
    end

    for con=1:length(electrode.contacts)
        vs=X*[electrode.contacts(con).vertices,ones(size(electrode.contacts(con).vertices,1),1)]';
        electrode.contacts(con).vertices=vs(1:3,:)';
        elrender=patch('Faces',electrode.contacts(con).faces,'Vertices',electrode.contacts(con).vertices);
        ea_specsurf(elrender,elspec.contact_color,aData);
    end

    axis equal
    view(0,0);
end
