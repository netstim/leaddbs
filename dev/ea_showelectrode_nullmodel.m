function electrode=ea_showelectrode_nullmodel(varargin)
% This function renders the electrode as defined by options.elspec and
% coords_mm.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if nargin
    options.elmodel=varargin{1};
else
    options.elmodel='Medtronic 3389';
end

pt=1;

options.sides=1;
elstruct.name=options.elmodel;
options=ea_resolve_elspec(options);
elspec=options.elspec;
resultfig=figure;

jetlist=othercolor('BuOr_12');
%   jetlist=jet;


for side=1:length(options.sides)
    %% nullmodel:
    coords_mm{side}=[0,0,1.5+0.75;0,0,1.5+0.75+1*2;0,0,1.5+0.75+2*2;0,0,1.5+0.75+3*2];
    trajectory{side}=[zeros(30,2),linspace(30,0,30)'];
    %%
    trajvector=mean(diff(trajectory{side}));
    trajvector=trajvector/norm(trajvector);


        startpoint=trajectory{side}(1,:);
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


    ea_specsurf(elrender{side}(1),usecolor,aData); ea_specsurf(elrender{side}(2),usecolor,aData); ea_specsurf(elrender{side}(3),usecolor,aData);

    cnt=4;

    % draw contacts
    for cntct=1:elspec.numel

        set(0,'CurrentFigure',resultfig);


        [elrender{side}(cnt),elrender{side}(cnt+1),elrender{side}(cnt+2)]=ea_cylinder(coords_mm{side}(cntct,:)-trajvector*(elspec.contact_length/2),coords_mm{side}(cntct,:)+trajvector*(elspec.contact_length/2),elspec.contact_diameter/2,100,repmat(elspec.contact_color,1,3),1,0);

        ea_specsurf(elrender{side}(cnt),elspec.contact_color,aData); ea_specsurf(elrender{side}(cnt+1),elspec.contact_color,aData); ea_specsurf(elrender{side}(cnt+2),elspec.contact_color,aData);
        cnt=cnt+3;
    end

    % draw trajectory between contacts
    for cntct=1:elspec.numel-1
        set(0,'CurrentFigure',resultfig);

        [elrender{side}(cnt),elrender{side}(cnt+1),elrender{side}(cnt+2)]=ea_cylinder(coords_mm{side}(cntct,:)-trajvector*(elspec.contact_length/2),coords_mm{side}(cntct+1,:)+trajvector*(elspec.contact_length/2),elspec.lead_diameter/2,100,repmat(elspec.lead_color,1,3),1,0);

        ea_specsurf(elrender{side}(cnt),usecolor,aData); ea_specsurf(elrender{side}(cnt+1),usecolor,aData); ea_specsurf(elrender{side}(cnt+2),usecolor,aData);
        cnt=cnt+3;
    end

    % draw tip

    if isfield(elstruct,'group')
        usecolor=elstruct.groupcolors(elstruct.group,:);
    else
        usecolor=elspec.tip_color;
    end
    set(0,'CurrentFigure',resultfig);

    [cX,cY,cZ] = cylinder((repmat(elspec.tip_diameter/2,1,10)-([10:-1:1].^10/10^10)*(elspec.tip_diameter/2)));

    cZ=cZ.*(elspec.tip_length); % scale to fit tip-diameter

    % define two points to define cylinder.
    X1=coords_mm{side}(1,:)+trajvector*(elspec.contact_length/2);
    X2=X1+trajvector*elspec.tip_length;


    cX=cX+X1(1);
    cY=cY+X1(2);
    cZ=cZ-(2*elspec.tip_length)/2+X1(3);



    elrender{side}(cnt)=surf(cX,cY,cZ);

    % Calulating the angle between the x direction and the required direction
    % of cylinder through dot product
    angle_X1X2=acos( dot( [0 0 -1],(X2-X1) )/( norm([0 0 -1])*norm(X2-X1)) )*180/pi;

    % Finding the axis of rotation (single rotation) to roate the cylinder in
    % X-direction to the required arbitrary direction through cross product
    axis_rot=cross([0 0 -1],(X2-X1) );

    if any(axis_rot) || angle_X1X2
        rotate(elrender{side}(cnt),axis_rot,angle_X1X2,X1)
    end
    ea_specsurf(elrender{side}(cnt),usecolor,aData);



end


if ~exist('elrender','var')
    elrender=nan;
end
axis equal
view(0,0);


%% build model spec:

cnt=1; cntcnt=1; inscnt=1;
ea_dispercent(0,'Exporting electrode components');

for comp=1:options.elspec.numel*2+1
    ea_dispercent(comp/(options.elspec.numel*2+1));
    try % shaft and contacts, here three surface components
        cyl=elrender{side}(cnt); top=elrender{side}(cnt+1); bottom=elrender{side}(cnt+2);
        cyl = surf2patch(cyl,'triangles');
        cyl.faces=[cyl.faces;top.Faces+length(cyl.vertices);bottom.Faces+length(cyl.vertices)+length(top.Vertices)];
        cyl.vertices=[cyl.vertices;top.Vertices;bottom.Vertices];
        cnt=cnt+3;
    catch % tip of the electrode, here only one surface component..
        cyl=elrender{side}(cnt);
        cyl = surf2patch(cyl,'triangles');
    end
    % this following method takes quite some time... even more importantly,
    % the info will be transfered from mesh to volume and lateron back to
    % mesh again. For now, this is still the most convenient method.
    if comp>1 && comp<options.elspec.numel+2 % these are the CONTACTS
        electrode.contacts(cntcnt)=cyl;
        cntcnt=cntcnt+1;
    else % these are the insulated shaft, tip and spacings..
        electrode.insulation(inscnt)=cyl;
        inscnt=inscnt+1;
    end
end

ea_dispercent(1,'end');

electrode.electrode_model=elstruct.name;
electrode.head_position=[0,0,options.elspec.tip_length+options.elspec.numel*options.elspec.contact_length+(options.elspec.numel-1)*options.elspec.contact_spacing];
electrode.numel=options.elspec.numel;
electrode.contact_color=options.elspec.contact_color;
electrode.lead_color=options.elspec.lead_color;
