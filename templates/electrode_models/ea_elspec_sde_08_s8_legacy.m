function electrode=ea_elspec_sde_08_s8_legacy(varargin)
% This function creates the electrode specification for a certain
% lead. Since this code is usually only executed once (to
% establish the model), it is not optimized in any way. You can however use
% this code to modify the electrode model and/or duplicate the function to
% build a different model.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

options.elmodel='SDE-08 S8 Legacy';

if nargin
    vizz=0;
else
    vizz=1;
end

pt=1;

options.sides=1;
elstruct.name=options.elmodel;
options=ea_resolve_elspec(options);
elspec=options.elspec;
resultfig=figure('visible','off');

jetlist=othercolor('BuOr_12');
%   jetlist=jet;
N=200; % resolution of electrode points

for iside=1:length(options.sides)
    side=options.sides(iside);
    %% nullmodel:
    coords_mm{side}=[0,0,elspec.tip_length+(elspec.contact_length/2);...
        0,0,elspec.tip_length+(elspec.contact_length/2)+1*(elspec.contact_spacing+elspec.contact_length);...
        0,0,elspec.tip_length+(elspec.contact_length/2)+2*(elspec.contact_spacing+elspec.contact_length);...
        0,0,elspec.tip_length+(elspec.contact_length/2)+3*(elspec.contact_spacing+elspec.contact_length);...
        0,0,elspec.tip_length+(elspec.contact_length/2)+4*(elspec.contact_spacing+elspec.contact_length);...
        0,0,elspec.tip_length+(elspec.contact_length/2)+5*(elspec.contact_spacing+elspec.contact_length);...
        0,0,elspec.tip_length+(elspec.contact_length/2)+6*(elspec.contact_spacing+elspec.contact_length);...
        0,0,elspec.tip_length+(elspec.contact_length/2)+7*(elspec.contact_spacing+elspec.contact_length)];

    trajectory{side}=[zeros(30,2),linspace(30,0,30)'];
    %%
    trajvector=mean(diff(trajectory{side}));
    trajvector=trajvector/norm(trajvector);


    startpoint=trajectory{side}(1,:)-(1.5*(coords_mm{side}(1,:)-trajectory{side}(1,:)));
    set(0,'CurrentFigure',resultfig);
    hold on
    % draw patientname
    lstartpoint=startpoint-(0.03*(coords_mm{side}(1,:)-startpoint));
    ellabel(side)=text(lstartpoint(1),lstartpoint(2),lstartpoint(3),elstruct.name);


    % draw trajectory
    lowerpoint=coords_mm{side}(elspec.numel,:)-trajvector*(elspec.contact_length/2);
    set(0,'CurrentFigure',resultfig);
    diams=repmat(elspec.lead_diameter/2,1,2);
    [cX,cY,cZ] = ea_singlecylinder((diams),N);

    cZ=cZ.*(startpoint(3)-lowerpoint(3)); % scale to fit tip-diameter
    cZ=cZ+lowerpoint(3);

    p=surf2patch(surf(cX,cY,cZ),'triangles');





    % add meshing-version to it
    [cX,cY,cZ] = ea_singlecylinder((diams),20);

    cZ=cZ.*(startpoint(3)-lowerpoint(3)); % scale to fit tip-diameter
    cZ=cZ+lowerpoint(3);
    a=surf2patch(surf(cX,cY,cZ));
    a=ea_reordercylinder(a,2);
    meshel.ins{1}.faces=a.faces;
    meshel.ins{1}.vertices=a.vertices;
    ndiv=length(meshel.ins{1}.vertices)/2;
    meshel.ins{1}.endplates=[1:ndiv;ndiv+1:2*ndiv];

    % add endplates:
    p.vertices=[p.vertices;...
        [0,0,startpoint(3)];...
        [0,0,lowerpoint(3)]];
    for pt=2:2:(N-1)*2
        p.faces=[p.faces;pt,pt+2,length(p.vertices)-1];
    end
    for pt=1:2:(N-1)*2
        p.faces=[p.faces;pt,pt+2,length(p.vertices)];
    end

    elrender{side}(1)=patch(p);

    cnt=2;



    if isfield(elstruct,'group')
        usecolor=elstruct.groupcolors(elstruct.group,:);
    else
        usecolor=elspec.lead_color;
    end


    aData=1;


    specsurf(elrender{side}(1),usecolor,aData);

    % draw contacts
    for cntct=1:elspec.numel
        set(0,'CurrentFigure',resultfig);
        diams=repmat(elspec.contact_diameter/2,1,2);
        [cX,cY,cZ] = ea_singlecylinder((diams),N);

        cZ=cZ.*(elspec.contact_length); % scale to fit tip-diameter
        htd=(max(cZ(:))/2);
        cZ=cZ-htd;
        cZ=cZ+coords_mm{side}(cntct,3);

        p=surf2patch(surf(cX,cY,cZ),'triangles');




        % add meshing-version to it
        [cX,cY,cZ] = ea_singlecylinder((diams),20);

        cZ=cZ.*(elspec.contact_length); % scale to fit tip-diameter
        htd=(max(cZ(:))/2);
        cZ=cZ-htd;
        cZ=cZ+coords_mm{side}(cntct,3);
        a=surf2patch(surf(cX,cY,cZ));

        a=ea_reordercylinder(a,2);
        meshel.con{cntct}.faces=a.faces;
        meshel.con{cntct}.vertices=a.vertices;
        ndiv=length(meshel.con{cntct}.vertices)/2;
        meshel.con{cntct}.endplates=[1:ndiv;ndiv+1:2*ndiv];


        % add endplates:
        p.vertices=[p.vertices;...
            coords_mm{side}(cntct,:)+[0,0,htd];...
            coords_mm{side}(cntct,:)-[0,0,htd]];
        for pt=2:2:(N-1)*2
            p.faces=[p.faces;pt,pt+2,length(p.vertices)-1];
        end
        for pt=1:2:(N-1)*2
            p.faces=[p.faces;pt,pt+2,length(p.vertices)];
        end

        elrender{side}(cnt)=patch(p);
        cnt=cnt+1;

    end

    % draw trajectory between contacts
    for cntct=1:elspec.numel-1
        set(0,'CurrentFigure',resultfig);
        diams=repmat(elspec.lead_diameter/2,1,2);
        [cX,cY,cZ] = ea_singlecylinder((diams),N);

        cZ=cZ.*(elspec.contact_spacing); % scale to fit tip-diameter
        htd=(max(cZ(:))/2);
        cZ=cZ-htd;
        hait=coords_mm{side}(cntct,3)+elspec.contact_length/2+elspec.contact_spacing/2;
        cZ=cZ+hait;

        p=surf2patch(surf(cX,cY,cZ),'triangles');


        % add meshing-version to it
        [cX,cY,cZ] = ea_singlecylinder((diams),20);

        cZ=cZ.*(elspec.contact_spacing); % scale to fit tip-diameter
        htd=(max(cZ(:))/2);
        cZ=cZ-htd;
        hait=coords_mm{side}(cntct,3)+elspec.contact_length/2+elspec.contact_spacing/2;
        cZ=cZ+hait;

        a=surf2patch(surf(cX,cY,cZ));
        a=ea_reordercylinder(a,2);
        meshel.ins{cntct+1}.faces=a.faces;
        meshel.ins{cntct+1}.vertices=a.vertices;
        ndiv=length(meshel.ins{cntct+1}.vertices)/2;
        meshel.ins{cntct+1}.endplates=[1:ndiv;ndiv+1:2*ndiv];


        % add endplates:
        p.vertices=[p.vertices;...
            [0,0,hait+htd];...
            [0,0,hait-htd]];
        for pt=2:2:(N-1)*2
            p.faces=[p.faces;pt,pt+2,length(p.vertices)-1];
        end
        for pt=1:2:(N-1)*2
            p.faces=[p.faces;pt,pt+2,length(p.vertices)];
        end

        elrender{side}(cnt)=patch(p);
        cnt=cnt+1;
    end








    % draw tip

    if isfield(elstruct,'group')
        usecolor=elstruct.groupcolors(elstruct.group,:);
    else
        usecolor=elspec.tip_color;
    end
    set(0,'CurrentFigure',resultfig);

    tipdiams=repmat(elspec.tip_diameter/2,1,19)-([10:-0.5:1].^10/10^10)*(elspec.tip_diameter/2);
    tipdiams(end+1)=elspec.tip_diameter/2;
    [cX,cY,cZ] = ea_singlecylinder((tipdiams),N);

    cZ=cZ.*(elspec.tip_length); % scale to fit tip-diameter

    % define two points to define cylinder.
    X1=coords_mm{side}(1,:)+trajvector*(elspec.contact_length/2);
    X2=X1+trajvector*elspec.tip_length;


    cX=cX+X1(1);
    cY=cY+X1(2);
    cZ=cZ-(2*elspec.tip_length)/2+X1(3);

    p=surf2patch(surf(cX,cY,cZ),'triangles');

    % add meshing-version to it

    [cX,cY,cZ] = ea_singlecylinder((tipdiams),20);

    cZ=cZ.*(elspec.tip_length); % scale to fit tip-diameter

    % define two points to define cylinder.
    X1=coords_mm{side}(1,:)+trajvector*(elspec.contact_length/2);
    X2=X1+trajvector*elspec.tip_length;


    cX=cX+X1(1);
    cY=cY+X1(2);
    cZ=cZ-(2*elspec.tip_length)/2+X1(3);
    a=surf2patch(surf(cX,cY,cZ));

    a=ea_reordercylinder(a,20);
    meshel.ins{end+1}.faces=a.faces;
    meshel.ins{end}.vertices=a.vertices;
    ndiv=100;
    meshel.ins{end}.endplates=[1:ndiv];


    % add endplate:
    p.vertices=[p.vertices;...
        [0,0,elspec.tip_length]];

    for pt=20:20:(length(p.vertices)-20)
        p.faces=[p.faces;pt,pt+20,length(p.vertices)];
    end



    elrender{side}(cnt)=patch(p);

    % Calulating the angle between the x direction and the required direction
    % of cylinder through dot product
    angle_X1X2=acos( dot( [0 0 -1],(X2-X1) )/( norm([0 0 -1])*norm(X2-X1)) )*180/pi;

    % Finding the axis of rotation (single rotation) to roate the cylinder in
    % X-direction to the required arbitrary direction through cross product
    axis_rot=cross([0 0 -1],(X2-X1) );

    if any(axis_rot) || angle_X1X2
        %       rotate(elrender{side}(cnt),axis_rot,angle_X1X2,X1)
    end
    specsurf(elrender{side}(cnt),usecolor,aData);



end




if ~exist('elrender','var')
    elrender=nan;
end
axis equal
view(0,0);




%% build model spec:

cnt=1; cntcnt=1; inscnt=1;
ea_dispercent(0,'Exporting electrode components');

for comp=1:elspec.numel*2+1
    ea_dispercent(comp/(elspec.numel*2+1));


    cyl=elrender{side}(cnt);
    cnt=cnt+1;


    if comp>1 && comp<elspec.numel+2 % these are the CONTACTS
        electrode.contacts(cntcnt).vertices=cyl.Vertices;
        electrode.contacts(cntcnt).faces=cyl.Faces;
        electrode.contacts(cntcnt).facevertexcdata=cyl.FaceVertexCData;
        cntcnt=cntcnt+1;
    else % these are the insulated shaft, tip and spacings..
        electrode.insulation(inscnt).vertices=cyl.Vertices;
        electrode.insulation(inscnt).faces=cyl.Faces;
        electrode.insulation(inscnt).facevertexcdata=cyl.FaceVertexCData;
        inscnt=inscnt+1;
    end
end

ea_dispercent(1,'end');

electrode.electrode_model=elstruct.name;
electrode.head_position=[0,0,elspec.tip_length+0.5*elspec.contact_length];
electrode.tail_position=[0,0,elspec.tip_length+4*elspec.contact_length+(4-1)*elspec.contact_spacing-0.5*elspec.contact_length];

electrode.x_position=[elspec.lead_diameter/2,0,elspec.tip_length+0.5*elspec.contact_length];
electrode.y_position=[0,elspec.lead_diameter/2,elspec.tip_length+0.5*elspec.contact_length];

electrode.numel=elspec.numel;
electrode.contact_color=elspec.contact_color;
electrode.lead_color=elspec.lead_color;
electrode.coords_mm=coords_mm{side};
electrode.meshel=meshel;
save([ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname],'electrode');

if vizz
    % visualize
    cnt=1;
    g=figure;
    X=eye(4);

    for ins=1:length(electrode.insulation)

        vs=X*[electrode.insulation(ins).vertices,ones(size(electrode.insulation(ins).vertices,1),1)]';
        electrode.insulation(ins).vertices=vs(1:3,:)';
        elrender{side}(cnt)=patch('Faces',electrode.insulation(ins).faces,'Vertices',electrode.insulation(ins).vertices);
        if isfield(elstruct,'group')
            usecolor=elstruct.groupcolors(elstruct.group,:);
        else
            usecolor=elspec.lead_color;
        end
        specsurf(elrender{side}(cnt),usecolor,aData);
        cnt=cnt+1;
    end
    for con=1:length(electrode.contacts)

        vs=X*[electrode.contacts(con).vertices,ones(size(electrode.contacts(con).vertices,1),1)]';
        electrode.contacts(con).vertices=vs(1:3,:)';
        elrender{side}(cnt)=patch('Faces',electrode.contacts(con).faces,'Vertices',electrode.contacts(con).vertices);

        specsurf(elrender{side}(cnt),elspec.contact_color,aData);

        cnt=cnt+1;
    end

    axis equal
    view(0,0);
end




%% build volumetric addition to it:
ea_genvol_sde_08_legacy(meshel,elspec,vizz);


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
%set(surfc,'FaceAlpha',0.3);

function C=rgb(C) % returns rgb values for the colors.

C = rem(floor((strfind('kbgcrmyw', C) - 1) * [0.25 0.5 1]), 2);
