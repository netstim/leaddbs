function electrode=ea_elspec_bostonvercisedirected(varargin)
% This function creates the electrode specification for a certain
% lead. Since this code is usually only executed once (to
% establish the model), it is not optimized in any way. You can however use
% this code to modify the electrode model and/or duplicate the function to
% build a different model.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


    options.elmodel='Boston Scientific Vercise Directed';

pt=1;

options.sides=1;
elstruct.name=options.elmodel;
options=ea_resolve_elspec(options);
elspec=options.elspec;
resultfig=figure('visible','off');

jetlist=othercolor('BuOr_12');
%   jetlist=jet;


for side=options.sides
    %% nullmodel:
    coords_mm{side}=[0,0,1.5+0.75;0,0,1.5+0.75+1*2;0,0,1.5+0.75+2*2;0,0,1.5+0.75+3*2];
    trajectory{side}=[zeros(30,2),linspace(30,0,30)'];
    %%
    trajvector=mean(diff(trajectory{side}));
    trajvector=trajvector/norm(trajvector);
    
    
    startpoint=trajectory{side}(1,:)-(1.5*(coords_mm{side}(1,:)-trajectory{side}(1,:)));
    set(0,'CurrentFigure',resultfig);
    
    % draw patientname
    lstartpoint=startpoint-(0.03*(coords_mm{side}(1,:)-startpoint));
    ellabel(side)=text(lstartpoint(1),lstartpoint(2),lstartpoint(3),elstruct.name);
    
    
    % draw trajectory
    [elrender{side}(1),elrender{side}(2),elrender{side}(3)]=ea_cylinder(startpoint,coords_mm{side}(elspec.numel,:)-trajvector*(elspec.contact_length/2),elspec.lead_diameter/2,100,repmat(elspec.lead_color,1,3),1,0);
    
    
    
    if isfield(elstruct,'group')
        usecolor=elstruct.groupcolors(elstruct.group,:);
    else
        usecolor=elspec.lead_color;
    end
    
    
    aData=1;
    
    
    specsurf(elrender{side}(1),usecolor,aData); specsurf(elrender{side}(2),usecolor,aData); specsurf(elrender{side}(3),usecolor,aData);
    
    cnt=4;
    
    % draw contacts
    for cntct=1:8 % first contact is the tip (see below).
        
        set(0,'CurrentFigure',resultfig);
        
        if cntct==1 % tip!
            
            % draw tip
            usecolor=elspec.tip_color;
            set(0,'CurrentFigure',resultfig);
            
            [cX,cY,cZ] = cylinder((repmat(elspec.tip_diameter/2,1,10)-([10:-1:1].^10/10^10)*(elspec.tip_diameter/2)));
            
            cZ=cZ.*(elspec.tip_length); % scale to fit tip-diameter
            
            % define two points to define cylinder.
            X1=coords_mm{side}(1,:)+trajvector*(elspec.tip_length/2);
            
            
            cX=cX+X1(1);
            cY=cY+X1(2);
            cZ=cZ+X1(3);
            
            elrender{side}(cnt)=surf(cX,cY,cZ);
            
            specsurf(elrender{side}(cnt),usecolor,aData);
            log(cnt)=1; % contact
            cnt=cnt+1;
        elseif cntct==2 || cntct==3
            
            scyl = ea_segmented_cylinder;
            for contact=1:3
                usecolor=elspec.contact_color;
                
                elrender{side}(cnt)=patch(scyl{1,contact});
                
                % scale size:
                elrender{side}(cnt).Vertices(:,1)=elrender{side}(cnt).Vertices(:,1).*(elspec.contact_diameter/2); % scale to fit tip-diameter
                elrender{side}(cnt).Vertices(:,2)=elrender{side}(cnt).Vertices(:,2).*(elspec.contact_diameter/2); % scale to fit tip-diameter
                elrender{side}(cnt).Vertices(:,3)=elrender{side}(cnt).Vertices(:,3).*(elspec.contact_length); % scale to fit tip-diameter
                
                % define two points to define cylinder.
                X1=coords_mm{side}(cntct,:)+trajvector*(elspec.contact_length/2);
                
                
                elrender{side}(cnt).Vertices(:,1)=elrender{side}(cnt).Vertices(:,1)+X1(1);
                elrender{side}(cnt).Vertices(:,2)=elrender{side}(cnt).Vertices(:,2)+X1(2);
                elrender{side}(cnt).Vertices(:,3)=elrender{side}(cnt).Vertices(:,3)+X1(3);
                specsurf(elrender{side}(cnt),usecolor,1);
                log(cnt)=1; % contact
                cnt=cnt+1;
            end
            
            
            
            for ins=1:4
                usecolor=elspec.lead_color;
                
                elrender{side}(cnt)=patch(scyl{2,ins});
                
                % scale size:
                elrender{side}(cnt).Vertices(:,1)=elrender{side}(cnt).Vertices(:,1).*(elspec.contact_diameter/2); % scale to fit tip-diameter
                elrender{side}(cnt).Vertices(:,2)=elrender{side}(cnt).Vertices(:,2).*(elspec.contact_diameter/2); % scale to fit tip-diameter
                elrender{side}(cnt).Vertices(:,3)=elrender{side}(cnt).Vertices(:,3).*(elspec.contact_length); % scale to fit tip-diameter
                
                % define two points to define cylinder.
                X1=coords_mm{side}(cntct,:)+trajvector*(elspec.contact_length/2);
                
                
                elrender{side}(cnt).Vertices(:,1)=elrender{side}(cnt).Vertices(:,1)+X1(1);
                elrender{side}(cnt).Vertices(:,2)=elrender{side}(cnt).Vertices(:,2)+X1(2);
                elrender{side}(cnt).Vertices(:,3)=elrender{side}(cnt).Vertices(:,3)+X1(3);
                specsurf(elrender{side}(cnt),usecolor,1);
                log(cnt)=0; % insulation
                
                cnt=cnt+1;
                
            end
            
            
            
        elseif cntct==4 % the only regular contact
            
            usecolor=elspec.lead_color;
            
            [elrender{side}(cnt),elrender{side}(cnt+1),elrender{side}(cnt+2)]=ea_cylinder(coords_mm{side}(cntct,:)-trajvector*(elspec.contact_length/2),coords_mm{side}(cntct,:)+trajvector*(elspec.contact_length/2),elspec.contact_diameter/2,100,repmat(elspec.contact_color,1,3),1,0);
            
            specsurf(elrender{side}(cnt),elspec.contact_color,aData); specsurf(elrender{side}(cnt+1),elspec.contact_color,aData); specsurf(elrender{side}(cnt+2),elspec.contact_color,aData);
            log(cnt:cnt+2)=1; % contact
            cnt=cnt+3;
            
        end
        
        
        
    end
    
    % draw trajectory between contacts
    for cntct=1:elspec.numel-1
        set(0,'CurrentFigure',resultfig);
        
        [elrender{side}(cnt),elrender{side}(cnt+1),elrender{side}(cnt+2)]=ea_cylinder(coords_mm{side}(cntct,:)-trajvector*(elspec.contact_length/2),coords_mm{side}(cntct+1,:)+trajvector*(elspec.contact_length/2),elspec.lead_diameter/2,100,repmat(elspec.lead_color,1,3),1,0);
        
        specsurf(elrender{side}(cnt),usecolor,aData); specsurf(elrender{side}(cnt+1),usecolor,aData); specsurf(elrender{side}(cnt+2),usecolor,aData);
        log(cnt:cnt+2)=0; % insulation
        
        cnt=cnt+3;
    end
    
    
    
    
    
    
    
    
    
    
    
    
end


if ~exist('elrender','var')
    elrender=nan;
end
axis equal
view(0,0);




%% build model spec:

cntcnt=1; inscnt=1;


%% order of components:
% 1-3: shaft (insulation)
% 4: tip (=contact1)
% 5: contact2
% 6: contact3
% 7: contact4
% 8-11: spacings
% 12: contact5
% 13: contact6
% 14: contact7
% 15-18: spacings
% 19-21: contact8
% 22-30: spacings

for comp=[1,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,22,25,28]
    if comp==1 % shaft
        cyl=elrender{side}(comp); top=elrender{side}(comp+1); bottom=elrender{side}(comp+2);
        cyl = surf2patch(cyl,'triangles');
        cyl.faces=[cyl.faces;top.Faces+length(cyl.vertices);bottom.Faces+length(cyl.vertices)+length(top.Vertices)];
        cyl.vertices=[cyl.vertices;top.Vertices;bottom.Vertices];
        electrode.insulation(inscnt)=cyl;
        inscnt=inscnt+1;
    elseif ismember(comp,4:7) % first 4 contacts
        cyl=elrender{side}(comp);
        try % if not already a patch..
            cyl = surf2patch(cyl,'triangles');
        end
        
        try
            electrode.contacts(cntcnt).faces=cyl.Faces;
            electrode.contacts(cntcnt).vertices=cyl.Vertices;
            electrode.contacts(cntcnt).facevertexcdata=cyl.FaceVertexCData;
        catch
            electrode.contacts(cntcnt).faces=cyl.faces;
            electrode.contacts(cntcnt).vertices=cyl.vertices;
            electrode.contacts(cntcnt).facevertexcdata=cyl.facevertexcdata;
        end
        cntcnt=cntcnt+1;
    elseif ismember(comp,8:11) % insulation
        cyl=elrender{side}(comp);
        try % if not already a patch..
            cyl = surf2patch(cyl,'triangles');
        end
        
        try
            electrode.insulation(inscnt).faces=cyl.Faces;
            electrode.insulation(inscnt).vertices=cyl.Vertices;
            electrode.insulation(inscnt).facevertexcdata=cyl.FaceVertexCData;
        catch
            electrode.insulation(inscnt).faces=cyl.faces;
            electrode.insulation(inscnt).vertices=cyl.vertices;
            electrode.insulation(inscnt).facevertexcdata=cyl.facevertexcdata;
        end
        
        inscnt=inscnt+1;
    elseif ismember(comp,12:14) % contacts
        cyl=elrender{side}(comp);
        try % if not already a patch..
            cyl = surf2patch(cyl,'triangles');
        end
        try
            electrode.contacts(cntcnt).faces=cyl.Faces;
            electrode.contacts(cntcnt).vertices=cyl.Vertices;
            electrode.contacts(cntcnt).facevertexcdata=cyl.FaceVertexCData;
        catch
            electrode.contacts(cntcnt).faces=cyl.faces;
            electrode.contacts(cntcnt).vertices=cyl.vertices;
            electrode.contacts(cntcnt).facevertexcdata=cyl.facevertexcdata;
        end
        cntcnt=cntcnt+1;
    elseif ismember(comp,15:18) % insulation
        cyl=elrender{side}(comp);
        try % if not already a patch..
            cyl = surf2patch(cyl,'triangles');
        end
        try
            electrode.insulation(inscnt).faces=cyl.Faces;
            electrode.insulation(inscnt).vertices=cyl.Vertices;
            electrode.insulation(inscnt).facevertexcdata=cyl.FaceVertexCData;
        catch
            electrode.insulation(inscnt).faces=cyl.faces;
            electrode.insulation(inscnt).vertices=cyl.vertices;
            electrode.insulation(inscnt).facevertexcdata=cyl.facevertexcdata;
        end
        inscnt=inscnt+1;
    elseif comp==19 % contact 8
        cyl=elrender{side}(comp); top=elrender{side}(comp+1); bottom=elrender{side}(comp+2);
        cyl = surf2patch(cyl,'triangles');
        cyl.faces=[cyl.faces;top.Faces+length(cyl.vertices);bottom.Faces+length(cyl.vertices)+length(top.Vertices)];
        cyl.vertices=[cyl.vertices;top.Vertices;bottom.Vertices];
        electrode.contacts(cntcnt)=cyl;
        cntcnt=cntcnt+1;
    elseif ismember(comp,[22,25,28]) % insulation
        cyl=elrender{side}(comp); top=elrender{side}(comp+1); bottom=elrender{side}(comp+2);
        cyl = surf2patch(cyl,'triangles');
        cyl.faces=[cyl.faces;top.Faces+length(cyl.vertices);bottom.Faces+length(cyl.vertices)+length(top.Vertices)];
        cyl.vertices=[cyl.vertices;top.Vertices;bottom.Vertices];
        electrode.insulation(inscnt)=cyl;
        inscnt=inscnt+1;
    end
    
end

electrode.electrode_model=elstruct.name;
electrode.head_position=[0,0,options.elspec.tip_length+0.5*options.elspec.contact_length];
electrode.tail_position=[0,0,options.elspec.tip_length+options.elspec.numel*options.elspec.contact_length+(options.elspec.numel-1)*options.elspec.contact_spacing-0.5*options.elspec.contact_length];

electrode.x_position=[options.elspec.lead_diameter/2,0,options.elspec.tip_length+0.5*options.elspec.contact_length];
electrode.y_position=[0,options.elspec.lead_diameter/2,options.elspec.tip_length+0.5*options.elspec.contact_length];

electrode.numel=8;
electrode.contact_color=options.elspec.contact_color;
electrode.lead_color=options.elspec.lead_color;

% add contact coordinates:
electrode.coords_mm(1,:)=coords_mm{side}(1,:);
electrode.coords_mm(2,:)=coords_mm{side}(2,:)+[-0.66,0,0];
electrode.coords_mm(3,:)=coords_mm{side}(2,:)+[0.33,0.66,0];
electrode.coords_mm(4,:)=coords_mm{side}(2,:)+[0.33,-0.66,0];
electrode.coords_mm(5,:)=coords_mm{side}(3,:)+[-0.66,0,0];
electrode.coords_mm(6,:)=coords_mm{side}(3,:)+[0.33,0.66,0];
electrode.coords_mm(7,:)=coords_mm{side}(3,:)+[0.33,-0.66,0];
electrode.coords_mm(8,:)=coords_mm{side}(4,:);

save([fileparts(which('lead')),filesep,'templates',filesep,'electrode_models',filesep,elspec.matfname],'electrode');
% visualize
cnt=1;
if ~nargin
g=figure;
else
    axes(varargin{1});
end
X=eye(4);

for ins=1:length(electrode.insulation)
    electrode.insulation(ins).vertices=X*[electrode.insulation(ins).vertices,ones(size(electrode.insulation(ins).vertices,1),1)]';
    electrode.insulation(ins).vertices=electrode.insulation(ins).vertices(1:3,:)';
    elrender{side}(cnt)=patch(electrode.insulation(ins));
    if isfield(elstruct,'group')
        usecolor=elstruct.groupcolors(elstruct.group,:);
    else
        usecolor=elspec.lead_color;
    end
    specsurf(elrender{side}(cnt),usecolor,aData);
    cnt=cnt+1;
end
for con=1:length(electrode.contacts)
    electrode.contacts(con).vertices=X*[electrode.contacts(con).vertices,ones(size(electrode.contacts(con).vertices,1),1)]';
    electrode.contacts(con).vertices=electrode.contacts(con).vertices(1:3,:)';
    elrender{side}(cnt)=patch(electrode.contacts(con));
    
    specsurf(elrender{side}(cnt),elspec.contact_color,aData);
    
    cnt=cnt+1;
end

axis equal
view(0,0);



function m=maxiso(cellinp) % simply returns the highest entry of matrices in a cell.
m=0;
for c=1:length(cellinp)
    nm=max(cellinp{c}(:));
    if nm>m; m=nm; end
end

function m=miniso(cellinp)
m=inf;
for c=1:length(cellinp)
    nm=min(cellinp{c}(:));
    if nm<m; m=nm; end
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
    vertices=get(surfc,'Vertices');
    cd=zeros(size(vertices));
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

function C=rgb(C) % returns rgb values for the colors.

C = rem(floor((strfind('kbgcrmyw', C) - 1) * [0.25 0.5 1]), 2);