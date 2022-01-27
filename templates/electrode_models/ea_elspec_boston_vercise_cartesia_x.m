function electrode = ea_elspec_boston_vercise_cartesia_x(vizz)
% Creates the FEM-compatible electrode model for Boston Scientific Vercise Cartesia X.
% Based on the mesh generated using Blender and tetgen, script created by Johannes Achtzehn.
% _________________________________________________________________________
% Copyright (C) 2022 Charite University Medicine Berlin
% Ningfei Li

% Set folder
elemodelPath = fileparts(mfilename('fullpath'));
modelFolder = 'Boston_Scientific_Vercise_Cartesia_X';

% Get specification
options.elmodel = 'Boston Scientific Vercise Cartesia X';
options = ea_resolve_elspec(options);
elspec = options.elspec;

% Get insulation and contact numbers
ins = ea_regexpdir([elemodelPath, filesep, modelFolder,filesep, 'Insulations'], '.*\.ele$');
numIns = numel(ins);
con = ea_regexpdir([elemodelPath, filesep, modelFolder,filesep, 'Contacts'], '.*\.ele$');
numCon = numel(con);

%% Import insulations and contacts meshes
for k = 1:numIns
    filename = ins{k}(1:end-4);
    [node,~,face] = readtetgen(filename);
    electrode.insulation(k).vertices = node;
    electrode.insulation(k).faces = face(:,1:3);
    clear face node filename
end

for k = 1:numCon
    filename = con{k}(1:end-4);
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
electrode.coords_mm = zeros(elspec.numel,3);

% Segmented contacts
for level=1:5
    levelHeight = elspec.tip_length*~elspec.tipiscontact + (level-1)*(elspec.contact_length+elspec.contact_spacing) + elspec.contact_length/2;
    levelIndStart = (level-1)*3 + 1;
    electrode.coords_mm(levelIndStart,:) = [0 0 levelHeight] + [0, elspec.lead_diameter/2, 0];
    electrode.coords_mm(levelIndStart+1,:) = [0 0 levelHeight] + [-cx, -cy, 0];
    electrode.coords_mm(levelIndStart+2,:) = [0 0 levelHeight] + [cx, -cy, 0];
end

% Normal contacts
for ind=levelIndStart+3:elspec.numel
    level = level + 1;
    electrode.coords_mm(ind,3) = elspec.tip_length*~elspec.tipiscontact + ...
                     elspec.contact_length/2 + ...
                     (level-1)*(elspec.contact_spacing+elspec.contact_length);
end

electrode.head_position = [0, 0, electrode.coords_mm(1,3)];
electrode.tail_position = [0, 0, electrode.coords_mm(10,3)];
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
    view(180,0);
end
