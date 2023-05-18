function electrode=ea_elspec_adtech_rd08r_sp05x(varargin)
% Generate electrode model of AdTech RD08R-SP05X

options.elmodel='AdTech RD08R-SP05X';

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

% rescale for tetgen model generation
rescaleratio = 50;
elspec.lead_diameter = elspec.lead_diameter*rescaleratio;
elspec.contact_length = elspec.contact_length*rescaleratio;
elspec.contact_diameter = elspec.contact_diameter*rescaleratio;
elspec.tip_diameter = elspec.tip_diameter*rescaleratio;
elspec.tip_length = elspec.tip_length*rescaleratio;
elspec.contact_spacing = elspec.contact_spacing*rescaleratio;

resultfig=figure('visible','off');
jetlist=othercolor('BuOr_12');
%   jetlist=jet;
N=10*rescaleratio; % resolution of electrode points

for iside=1:length(options.sides)
    side=options.sides(iside);
    % nullmodel:
    coords_mm{side}=[0,0,elspec.tip_length+(elspec.contact_length/2);...
        0,0,elspec.tip_length+(elspec.contact_length/2)+1*(elspec.contact_spacing+elspec.contact_length);...
        0,0,elspec.tip_length+(elspec.contact_length/2)+2*(elspec.contact_spacing+elspec.contact_length);...
        0,0,elspec.tip_length+(elspec.contact_length/2)+3*(elspec.contact_spacing+elspec.contact_length);...
        0,0,elspec.tip_length+(elspec.contact_length/2)+4*(elspec.contact_spacing+elspec.contact_length);...
        0,0,elspec.tip_length+(elspec.contact_length/2)+5*(elspec.contact_spacing+elspec.contact_length);...
        0,0,elspec.tip_length+(elspec.contact_length/2)+6*(elspec.contact_spacing+elspec.contact_length);...
        0,0,elspec.tip_length+(elspec.contact_length/2)+7*(elspec.contact_spacing+elspec.contact_length);];

    trajectory{side}=[zeros(50*rescaleratio,2),linspace(50*rescaleratio,0,50*rescaleratio)'];
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

    ea_specsurf(elrender{side}(1),usecolor,aData);

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
    ea_specsurf(elrender{side}(cnt),usecolor,aData);
end

if ~exist('elrender','var')
    elrender=nan;
end
axis equal
view(0,0);

% build model spec:

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
        ea_specsurf(elrender{side}(cnt),usecolor,aData);
        cnt=cnt+1;
    end
    for con=1:length(electrode.contacts)

        vs=X*[electrode.contacts(con).vertices,ones(size(electrode.contacts(con).vertices,1),1)]';
        electrode.contacts(con).vertices=vs(1:3,:)';
        elrender{side}(cnt)=patch('Faces',electrode.contacts(con).faces,'Vertices',electrode.contacts(con).vertices);

        ea_specsurf(elrender{side}(cnt),elspec.contact_color,aData);

        cnt=cnt+1;
    end

    axis equal
    view(0,0);
end

%% build volumetric addition to it:
ea_genvol_adtech(meshel,elspec,rescaleratio,vizz);

%% Revert to real dimension
electrode.head_position = electrode.head_position/rescaleratio;
electrode.tail_position = electrode.tail_position/rescaleratio;
electrode.x_position = electrode.x_position/rescaleratio;
electrode.y_position = electrode.y_position/rescaleratio;
electrode.coords_mm = electrode.coords_mm/rescaleratio;

for i=1:length(electrode.insulation)
    electrode.insulation(i).vertices = electrode.insulation(i).vertices/rescaleratio;
end

for i=1:length(electrode.contacts)
    electrode.contacts(i).vertices = electrode.contacts(i).vertices/rescaleratio;
end

for i=1:length(meshel.ins)
    meshel.ins{i}.vertices = meshel.ins{1}.vertices * rescaleratio;
end

for i=1:length(meshel.con)
    meshel.con{i}.vertices = meshel.con{1}.vertices * rescaleratio;
end

save([ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname],'electrode');
