function electrode=ea_elspec_elaine_rat_electrode(varargin)
% This function creates the electrode specification for a certain
% lead. Since this code is usually only executed once (to
% establish the model), it is not optimized in any way. You can however use
% this code to modify the electrode model and/or duplicate the function to
% build a different model.

elemodelPath = fileparts(mfilename('fullpath'));

options.elmodel='ELAINE Rat Electrode';

%% import insulations and contacts from subfolder
for k = 1
    filename = [elemodelPath, filesep, 'ELAINE_Rat_Electrode_Components', ...
    	filesep, 'Insulations', filesep, 'ins', num2str(k), '.1'];
    [node,~,face]=readtetgen(filename);
    electrode.insulation(k).vertices = node;
    electrode.insulation(k).faces = face(:,1:3);
    clear face node filename
end

for k = 1
    filename = [elemodelPath, filesep, 'ELAINE_Rat_Electrode_Components', ...
    	filesep, 'Contacts', filesep, 'con', num2str(k), '.1'];
    [node,~,face]=readtetgen(filename);
    electrode.contacts(k).vertices = node;
    electrode.contacts(k).faces = face(:,1:3);
    clear face node filename
end

%% other specifications
options = ea_resolve_elspec(options);
elspec = options.elspec;

electrode.electrode_model = options.elmodel;
electrode.head_position = [0 0 elspec.tip_length/2]; % dummy value
electrode.tail_position = [0 0 elspec.tip_length/2+elspec.tip_length]; % dummy value
electrode.x_position = [elspec.lead_diameter/2, 0, elspec.tip_length/2]; % dummy value
electrode.y_position = [0, elspec.lead_diameter/2, elspec.tip_length/2]; % dummy value
electrode.numel = 0;
electrode.contact_color = 0.3;
electrode.lead_color = 0.7;

electrode.coords_mm(1,:)=[0 0 elspec.tip_length/2];

electrode.isdirected = 0;

%% saving electrode struct
save([elemodelPath, filesep, 'elaine_rat_electrode.mat'],'electrode');

%% create and save _vol file
filename = [elemodelPath, filesep, 'ELAINE_Rat_Electrode_Components', filesep, 'final.1'];
[node,~,face] = readtetgen(filename);
save([elemodelPath, filesep, 'elaine_rat_electrode_vol.mat'],'face','node')
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
