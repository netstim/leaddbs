function electrode=ea_elspec_aleva_directstim_directed(varargin)
% Creates the FEM-compatible electrode model for Aleva directSTIM.
% It's based on the mesh generated using SketchUp and tetgen.

% Set folder
elemodelPath = fileparts(mfilename('fullpath'));
modelFolder = 'Aleva_directSTIM_Directed';

% Get specification
options.elmodel = 'Aleva directSTIM Directed';
options = ea_resolve_elspec(options);     
elspec = options.elspec;

% Get insulation and contact numbers
numIns = numel(ea_regexpdir([elemodelPath, filesep, modelFolder,filesep, 'Insulations'], '.*\.face$'));
numCon = numel(ea_regexpdir([elemodelPath, filesep, modelFolder,filesep, 'Contacts'], '.*\.face$'));

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

% The segmented contacts are clockwise arranged seen from the top
% view, the same as in the models in the components folder.
electrode.coords_mm(1,:) = [0 0 1.82]+[0, elspec.lead_diameter/2, 0];
electrode.coords_mm(2,:) = [0 0 1.82]+[cx, -cy, 0];
electrode.coords_mm(3,:) = [0 0 1.82]+[-cx, -cy, 0];
electrode.coords_mm(4,:) = [0 0 3.82]+[0, elspec.lead_diameter/2, 0];
electrode.coords_mm(5,:) = [0 0 3.82]+[cx, -cy, 0];
electrode.coords_mm(6,:) = [0 0 3.82]+[-cx,-cy, 0];
electrode.coords_mm(7,:) = [0 0 5.82]+[0, elspec.lead_diameter/2, 0];
electrode.coords_mm(8,:) = [0 0 5.82]+[cx,-cy, 0];
electrode.coords_mm(9,:) = [0 0 5.82]+[-cx,-cy, 0];
electrode.coords_mm(10,:) = [0 0 7.82]+[0, elspec.lead_diameter/2, 0];
electrode.coords_mm(11,:) = [0 0 7.82]+[cx,-cy, 0];
electrode.coords_mm(12,:) = [0 0 7.82]+[-cx,-cy, 0];

electrode.head_position = [0 0 1.82];
electrode.tail_position = [0 0 7.82];        
electrode.x_position = [elspec.lead_diameter/2, 0, 1.82]; 
electrode.y_position = [0, elspec.lead_diameter/2, 1.82];

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
