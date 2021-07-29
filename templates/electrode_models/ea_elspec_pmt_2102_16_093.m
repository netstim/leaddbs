function electrode = ea_elspec_pmt_2102_16_093(vizz)
% Creates the FEM-compatible electrode model for PMT 2102-16-093.
% It's based on the mesh generated using SketchUp and tetgen.
% _________________________________________________________________________
% Copyright (C) 2020 Charite University Medicine Berlin
% Ningfei Li

% Set folder
elemodelPath = fileparts(mfilename('fullpath'));
modelFolder = 'PMT_2102-16-093';

% Get specification
options.elmodel = 'PMT 2102-16-093';
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
electrode.coords_mm = zeros(numCon,3);
for i=1:numCon
    electrode.coords_mm(i,3) = elspec.tip_length*~elspec.tipiscontact + ...
                     elspec.contact_length/2 + ...
                     (i-1)*(elspec.contact_spacing+elspec.contact_length);
end

electrode.head_position = electrode.coords_mm(1,:);
electrode.tail_position = electrode.coords_mm(4,:);
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
