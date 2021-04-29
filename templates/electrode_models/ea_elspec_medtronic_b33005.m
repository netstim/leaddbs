function electrode=ea_elspec_medtronic_b33005(varargin)
% This function creates the electrode specification for the Vercise directed lead.
% In contrast to other electrode generation functions, it is based on an
% externally generated FEM-compatible model, stored in the
% Medtronic_B33005 subfolder.
% __________________________________________________________________________________
% Copyright (C) 2021 Charite University Medicine Berlin, Movement Disorders Unit
% Ningfei Li

% Set folder
elemodelPath = fileparts(mfilename('fullpath'));
modelFolder = 'Medtronic_B33005';

% Get specification
options.elmodel = 'Medtronic B33005';
options = ea_resolve_elspec(options);
elspec = options.elspec;

% Get insulation and contact numbers
numIns = numel(ea_regexpdir([elemodelPath, filesep, modelFolder,filesep, 'Insulations'], '.*\.smesh$'));
numCon = numel(ea_regexpdir([elemodelPath, filesep, modelFolder,filesep, 'Contacts'], '.*\.smesh$'));

%% Import insulations and contacts meshes
for k = 1:numIns
    filename = [elemodelPath, filesep, modelFolder, filesep, 'Insulations', filesep, 'ins', num2str(k), '.1'];
    [node,~,face] = readtetgen(filename);
    electrode.insulation(k).vertices = node;
    electrode.insulation(k).faces = face(:,1:3);
    clear face node filename
end

for k = 1:numCon
    filename = [elemodelPath, filesep, modelFolder, filesep, 'Contacts', filesep, 'con', num2str(k), '.1'];
    [node,~,face] = readtetgen(filename);
    electrode.contacts(k).vertices = node;
    electrode.contacts(k).faces = face(:,1:3);
    clear face node filename
end

%% Contact coordinates and other specifications
cx = elspec.lead_diameter/2*cos(pi/6);
cy = elspec.lead_diameter/2*sin(pi/6);

% The segmented contacts are counter-clockwise arranged seen from the top
% view, the same as in the models in the components folder.
electrode.coords_mm(1,:) = [0 0 1.65];
electrode.coords_mm(2,:) = [0 0 3.65]+[0, elspec.lead_diameter/2, 0];
electrode.coords_mm(3,:) = [0 0 3.65]+[-cx, -cy, 0];
electrode.coords_mm(4,:) = [0 0 3.65]+[cx, -cy, 0];
electrode.coords_mm(5,:) = [0 0 5.65]+[0, elspec.lead_diameter/2, 0];
electrode.coords_mm(6,:) = [0 0 5.65]+[-cx,-cy, 0];
electrode.coords_mm(7,:) = [0 0 5.65]+[cx,-cy, 0];
electrode.coords_mm(8,:) = [0 0 7.75];

electrode.head_position = electrode.coords_mm(1,:);
electrode.tail_position = electrode.coords_mm(end,:);
electrode.x_position = [elspec.lead_diameter/2, 0, electrode.coords_mm(1,3)];
electrode.y_position = [0, elspec.lead_diameter/2, electrode.coords_mm(1,3)];

electrode.electrode_model = options.elmodel;
electrode.isdirected = elspec.isdirected;
electrode.numel = elspec.numel;
electrode.contact_color = elspec.contact_color;
electrode.lead_color = elspec.lead_color;

%% Save electrode model
save([elemodelPath, filesep, elspec.matfname, '.mat'], 'electrode');

%% Create and save *_vol.mat
filename = [elemodelPath, filesep, modelFolder, filesep, 'final.1'];
[node,~,face] = readtetgen(filename);
save([elemodelPath, filesep, elspec.matfname, '_vol.mat'], 'face', 'node')
clear node face

%% Visualize
if ~exist('vizz', 'var')
    vizz = 1;
end

if vizz
    X = eye(4);
    aData = 1;

    figure;
    for ins=1:length(electrode.insulation)
        vs = X*[electrode.insulation(ins).vertices,ones(size(electrode.insulation(ins).vertices,1),1)]';
        electrode.insulation(ins).vertices = vs(1:3,:)';
        elrender = patch('Faces',electrode.insulation(ins).faces,'Vertices',electrode.insulation(ins).vertices);
        ea_specsurf(elrender,elspec.lead_color,aData);
    end

    for con=1:length(electrode.contacts)
        vs = X*[electrode.contacts(con).vertices,ones(size(electrode.contacts(con).vertices,1),1)]';
        electrode.contacts(con).vertices = vs(1:3,:)';
        elrender = patch('Faces',electrode.contacts(con).faces,'Vertices',electrode.contacts(con).vertices);
        ea_specsurf(elrender,elspec.contact_color,aData);
    end

    axis equal
    view(0,0);
end
