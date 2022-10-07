function electrode=ea_elspec_abbott_directed_05(varargin)

% This function creates the electrode specification for the short Abbott directed lead.
% In contrast to other electrode generation functions, it is based on an
% externally generated FEM-compatible model, stored in the
% Abbott_Directed_05 subfolder.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

% Set folder
elemodelPath = fileparts(mfilename('fullpath'));
modelFolder = 'Abbott_Directed_05';

% Get insulation and contact numbers
numIns = numel(ea_regexpdir([elemodelPath, filesep, modelFolder,filesep, 'Insulations'], '.*\.smesh$'));
numCon = numel(ea_regexpdir([elemodelPath, filesep, modelFolder,filesep, 'Contacts'], '.*\.smesh$'));

%% import insulations and contacts from subfolder
for k = 1:numIns
    filename = [elemodelPath,filesep,modelFolder,filesep,'Insulations',filesep,'ins',num2str(k),'.1'];
    [node,~,face]=readtetgen(filename);
    electrode.insulation(k).vertices = node;
    electrode.insulation(k).faces = face(:,1:3);
    clear face node filename
end

for k = 1:numCon
    filename = [elemodelPath,filesep,modelFolder,filesep,'Contacts',filesep,'con',num2str(k),'.1'];
    [node,~,face]=readtetgen(filename);
    electrode.contacts(k).vertices = node;
    electrode.contacts(k).faces = face(:,1:3);
    clear face node filename
end

%% other specifications
options.elmodel = 'Abbott Directed 6172 (short)';
options = ea_resolve_elspec(options);
elspec = options.elspec;

electrode.electrode_model = options.elmodel;
electrode.head_position = [0 0 1.75];
electrode.tail_position = [0 0 7.75];
electrode.x_position = [elspec.lead_diameter/2 0 1.75];
electrode.y_position = [0 elspec.lead_diameter/2 1.75];
electrode.numel = 8;
electrode.contact_color = 0.3;
electrode.lead_color = 0.7;

cx = elspec.lead_diameter/2*cos(pi/6);
cy = elspec.lead_diameter/2*sin(pi/6);

electrode.coords_mm(1,:)=[0 0 1.75];
electrode.coords_mm(2,:)=[0 0 3.75]+[0, elspec.lead_diameter/2, 0];
electrode.coords_mm(3,:)=[0 0 3.75]+[cx, -cy, 0];
electrode.coords_mm(4,:)=[0 0 3.75]+[-cx, -cy, 0];
electrode.coords_mm(5,:)=[0 0 5.75]+[0, elspec.lead_diameter/2, 0];
electrode.coords_mm(6,:)=[0 0 5.75]+[cx, -cy, 0];
electrode.coords_mm(7,:)=[0 0 5.75]+[-cx, -cy, 0];
electrode.coords_mm(8,:)=[0 0 7.75];

electrode.isdirected = 1;

%% saving electrode struct
save([elemodelPath,filesep,'abbott_directed_05.mat'],'electrode');

%% create and save _vol file
filename = [elemodelPath,filesep,modelFolder,filesep,'final.1'];
[node,~,face]=readtetgen(filename);
save([elemodelPath,filesep,'abbott_directed_05_vol.mat'],'face','node')
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
