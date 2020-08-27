function electrode=ea_elspec_boston_vercise_directed(varargin)

% This function creates the electrode specification for the Vercise directed lead.
% In contrast to other electrode generation functions, it is based on an
% externally generated FEM-compatible model, stored in the
% Boston_Vercise_Directed_Components subfolder.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

elemodelPath = fileparts(mfilename('fullpath'));

% The segmented contacts are anti-clockwise arranged seen from the top view.
% But they are clockwise ordered in the models in the components folder. So
% We need to reoder it here.
electrodeorder = [1 2 4 3 5 7 6 8 9];

%% import insulations and contacts from subfolder
for k = 1:16
    filename = [elemodelPath, filesep, 'Boston_Vercise_Directed_Components', ...
    	filesep, 'Insulations', filesep, 'ins', num2str(k), '.1'];
    [node,~,face]=readtetgen(filename);
    node(:,1) = -node(:,1); % Flip X axis
    electrode.insulation(k).vertices = node;
    electrode.insulation(k).faces = face(:,1:3);
    clear face node filename
end

for k = 1:numel(electrodeorder)
    filename = [elemodelPath, filesep, 'Boston_Vercise_Directed_Components', ...
    	filesep, 'Contacts', filesep, 'con', num2str(electrodeorder(k)), '.1'];
    [node,~,face]=readtetgen(filename);
    node(:,1) = -node(:,1); % Flip X axis
    electrode.contacts(k).vertices = node;
    electrode.contacts(k).faces = face(:,1:3);
    clear face node filename
end

%% other specifications
electrode.electrode_model = 'Boston Scientific Vercise Directed';
electrode.head_position = [0 0 0.75];
electrode.tail_position = [0 0 6.75];
electrode.x_position = [0.65 0 0.75];
electrode.y_position = [0 0.65 0.75];
electrode.numel = 8;
electrode.contact_color = 0.3;
electrode.lead_color = 0.7;

% The segmented contact in the null model are also anti-clockwise arranged
% from the top view.
electrode.coords_mm(1,:)=[0 0 0.75];
electrode.coords_mm(2,:)=[0 0 2.75]+[-0.66,0,0];
electrode.coords_mm(3,:)=[0 0 2.75]+[0.33,-0.66,0];
electrode.coords_mm(4,:)=[0 0 2.75]+[0.33,0.66,0];
electrode.coords_mm(5,:)=[0 0 4.75]+[-0.66,0,0];
electrode.coords_mm(6,:)=[0 0 4.75]+[0.33,-0.66,0];
electrode.coords_mm(7,:)=[0 0 4.75]+[0.33,0.66,0];
electrode.coords_mm(8,:)=[0 0 6.75];

%% saving electrode struct
save([elemodelPath, filesep, 'boston_vercise_directed.mat'],'electrode');

%% create and save _vol file
filename = [elemodelPath, filesep, 'Boston_Vercise_Directed_Components', filesep, 'final.1'];
[node,~,face] = readtetgen(filename);
save([elemodelPath, filesep, 'boston_vercise_directed_vol.mat'],'face','node')
clear node face

%% visualize
if nargin
    vizz=0;
else
    vizz=1;
end

options.elmodel = 'Boston Scientific Vercise Directed';
options = ea_resolve_elspec(options);
elspec = options.elspec;

if vizz
    X = eye(4);
    aData = 1;

    figure;
    for ins=1:length(electrode.insulation)
        vs=X*[electrode.insulation(ins).vertices,ones(size(electrode.insulation(ins).vertices,1),1)]';
        electrode.insulation(ins).vertices=vs(1:3,:)';
        elrender=patch('Faces',electrode.insulation(ins).faces,'Vertices',electrode.insulation(ins).vertices);
        specsurf(elrender,elspec.lead_color,aData);
    end

    for con=1:length(electrode.contacts)
        vs=X*[electrode.contacts(con).vertices,ones(size(electrode.contacts(con).vertices,1),1)]';
        electrode.contacts(con).vertices=vs(1:3,:)';
        elrender=patch('Faces',electrode.contacts(con).faces,'Vertices',electrode.contacts(con).vertices);
        specsurf(elrender,elspec.contact_color,aData);
    end

    axis equal
    view(0,0);
end


function specsurf(varargin)

surfc=varargin{1};
color=varargin{2};
if nargin==3
    aData=varargin{3};
end

len=get(surfc,'ZData');

cd=zeros([size(len),3]);
cd(:,:,1)=color(1);
try % works if color is denoted as 1x3 array
    cd(:,:,2)=color(2);cd(:,:,3)=color(3);
catch % if color is denoted as gray value (1x1) only
    cd(:,:,2)=color(1);cd(:,:,3)=color(1);
end

cd=cd+0.01*randn(size(cd));

set(surfc,'FaceColor','interp');
set(surfc,'CData',cd);

try % for patches
    Vertices=get(surfc,'Vertices');
    cd=zeros(size(Vertices));
    cd(:)=color(1);
    set(surfc,'FaceVertexCData',cd);
end
set(surfc,'AlphaDataMapping','none');

set(surfc,'FaceLighting','phong');
set(surfc,'SpecularColorReflectance',0);
set(surfc,'SpecularExponent',10);
set(surfc,'EdgeColor','none')

if nargin==3
    set(surfc,'FaceAlpha',aData);
end
