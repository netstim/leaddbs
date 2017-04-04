function electrode=ea_elspec_boston_vercise_directed(varargin)
% This function creates the electrode specification for a certain
% lead. Since this code is usually only executed once (to
% establish the model), it is not optimized in any way. You can however use
% this code to modify the electrode model and/or duplicate the function to
% build a different model.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if nargin
    vizz=0;
else
    vizz=1;
end
options.elmodel='Boston Scientific Vercise Directed';

pt=1;

options.sides=1;
elstruct.name=options.elmodel;
options=ea_resolve_elspec(options);
elspec=options.elspec;
resultfig=figure('visible','off');

jetlist=jet;
%   jetlist=jet;
N=200; % resolution of electrode points
aData=1;

cnt=4;
for side=1:length(options.sides)
    %% nullmodel:
    coords_mm{side}=[0,0,0+0.75;0,0,0+0.75+1*2;0,0,0+0.75+2*2;0,0,0+0.75+3*2];

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
    lowerpoint=coords_mm{side}(4,:)-trajvector*(elspec.contact_length/2);
    set(0,'CurrentFigure',resultfig);
    diams=repmat(elspec.lead_diameter/2,1,2);
    [cX,cY,cZ] = ea_singlecylinder((diams),N);

    cZ=cZ.*(startpoint(3)-lowerpoint(3)); % scale to fit tip-diameter
    cZ=cZ+lowerpoint(3);

    p=surf2patch(surf(cX,cY,cZ),'triangles');


    % add meshing-version to shaft
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



    % draw contacts
    for cntct=1:8 % first contact is the tip (see below).

        set(0,'CurrentFigure',resultfig);

        if cntct==1 % tip!

            % draw tip


            usecolor=elspec.tip_color;
            set(0,'CurrentFigure',resultfig);

            tipdiams=repmat(elspec.tip_diameter/2,1,19)-([10:-0.5:1].^10/10^10)*(elspec.tip_diameter/2);
            tipdiams(end+1)=elspec.tip_diameter/2;
            [cX,cY,cZ] = ea_singlecylinder((tipdiams),N);

            cZ=cZ.*(elspec.tip_length); % scale to fit tip-diameter

            % define two points to define cylinder.
            X2=coords_mm{side}(1,:)+trajvector*(elspec.contact_length/2);
            X1=X2-trajvector*elspec.tip_length;


            cX=cX+X1(1);
            cY=cY+X1(2);
            cZ=cZ-(2*elspec.tip_length)/2+X1(3);

            p=surf2patch(surf(cX,cY,cZ),'triangles');

            % add meshing-version to it

            [cX,cY,cZ] = ea_singlecylinder((tipdiams),20);

            cZ=cZ.*(elspec.tip_length); % scale to fit tip-diameter



            cX=cX+X1(1);
            cY=cY+X1(2);
            cZ=cZ-(2*elspec.tip_length)/2+X1(3);
            a=surf2patch(surf(cX,cY,cZ));

            a=ea_reordercylinder(a,20);
            meshel.con{1}.faces=a.faces;
            meshel.con{1}.vertices=a.vertices;
            ndiv=100;
            meshel.con{1}.endplates=[1:ndiv];


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


        elseif cntct==2 || cntct==3
            % define two points to define cylinder.
            X1=coords_mm{side}(cntct,:)+trajvector*(elspec.contact_length/2);
            [no,fc,seeds] = ea_segmented_cylinder_qq(60,1,0.5,1,3,0.8);
            [vnode,velem,vface]=s2m(no,fc,1,3);


            vcnt=1;
            for contact=1:3
                usecolor=elspec.contact_color;
                thisface=volface(velem(velem(:,end)==vcnt,1:4));
                [thisno,thisfc]=removeisolatednode(vnode,thisface);
                elrender{side}(cnt)=patch('Vertices',thisno,'Faces',thisfc);


                % scale size:
                elrender{side}(cnt).Vertices(:,1)=elrender{side}(cnt).Vertices(:,1).*(elspec.contact_diameter/2); % scale to fit tip-diameter
                elrender{side}(cnt).Vertices(:,2)=elrender{side}(cnt).Vertices(:,2).*(elspec.contact_diameter/2); % scale to fit tip-diameter
                elrender{side}(cnt).Vertices(:,3)=elrender{side}(cnt).Vertices(:,3).*(elspec.contact_length); % scale to fit tip-diameter



                elrender{side}(cnt).Vertices(:,1)=elrender{side}(cnt).Vertices(:,1)+X1(1);
                elrender{side}(cnt).Vertices(:,2)=elrender{side}(cnt).Vertices(:,2)+X1(2);
                elrender{side}(cnt).Vertices(:,3)=elrender{side}(cnt).Vertices(:,3)+X1(3);
                specsurf(elrender{side}(cnt),usecolor,1);

                % meshing version, here just a duplicate
                meshel.con{end+1}.faces=elrender{side}(cnt).Faces;
                meshel.con{end}.vertices=elrender{side}(cnt).Vertices;
                if vizz
                    %                     figure
                    %                     hold on
                    %                         plotmesh(meshel.con{end}.vertices,meshel.con{end}.faces)
                end

                log(cnt)=1; % contact
                cnt=cnt+1;
                vcnt=vcnt+1;
            end



            for ins=1:4
                usecolor=elspec.lead_color;

                thisface=volface(velem(velem(:,end)==vcnt,1:4));
                [thisno,thisfc]=removeisolatednode(vnode,thisface);
                elrender{side}(cnt)=patch('Vertices',thisno,'Faces',thisfc);
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

                % meshing version, here just a duplicate
                    meshel.ins{end+1}.faces=elrender{side}(cnt).Faces;
                    meshel.ins{end}.vertices=elrender{side}(cnt).Vertices;



                log(cnt)=0; % insulation

                cnt=cnt+1;
                vcnt=vcnt+1;

            end



            % plotmesh(contno,contfc);
            % subplot(122);
            % plotmesh(insuno,insufc);
            %
            % figure
            % plotmesh(vnode,vface)






        elseif cntct==4 % the only regular contact

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
            meshel.con{end+1}.faces=a.faces;
            meshel.con{end}.vertices=a.vertices;
            ndiv=length(meshel.con{end}.vertices)/2;
            meshel.con{end}.endplates=[1:ndiv;ndiv+1:2*ndiv];


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



    end


    % draw trajectory between contacts

    for cntct=1:3

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

        meshel.ins{end+1}.faces=a.faces;
        meshel.ins{end}.vertices=a.vertices;
        ndiv=length(meshel.ins{end}.vertices)/2;
        meshel.ins{end}.endplates=[1:ndiv;ndiv+1:2*ndiv];


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






end


if ~exist('elrender','var')
    elrender=nan;
end
axis equal
view(0,0);




%% build model spec:

cntcnt=1; inscnt=1;



if vizz
    figure
    hold on

    for c=1:length(meshel.con)
        plotmesh(meshel.con{c}.vertices,meshel.con{c}.faces)
    end
    for c=1:length(meshel.ins)
        plotmesh(meshel.ins{c}.vertices,meshel.ins{c}.faces)
    end
end

% export vol version:

ea_genvol_bscdirected(meshel,elspec,1);

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
        try
            cyl = surf2patch(cyl,'triangles');
        catch
            keyboard
        end
        cyl.faces=[cyl.faces;top.Faces+length(cyl.vertices);bottom.Faces+length(cyl.vertices)+length(top.Vertices)];
        cyl.vertices=[cyl.vertices;top.Vertices;bottom.Vertices];
        electrode.insulation(inscnt)=cyl;
        inscnt=inscnt+1;



    elseif ismember(comp,4:7) % first 4 contacts
        cyl=elrender{side}(comp);

        try % if not already a patch..
            cyl = surf2patch(cyl,'triangles');
        end

        if comp==4 % tip
            % close lid (add endplate)


            % find all vertices at Z=3mm
            z3=cyl.vertices(cyl.vertices(:,3)==3,:);
            z3ix=find(cyl.vertices(:,3)==3);

            % add 0,0,3 (lid midpoint)
            cyl.vertices=[cyl.vertices;0,0,3];

            midpointix=length(cyl.vertices);

            for face=1:length(z3ix)
                try
                    cyl.faces=[cyl.faces;z3ix(face),z3ix(face+1),midpointix];
                catch
                    cyl.faces=[cyl.faces;z3ix(face),z3ix(1),midpointix];
                end
            end

            %
            %
            %             z3ix=z3ix';
            %             keyboard
            %             vtc=nan(length(cyl.faces)+1,21);
            %             vtc(1:length(cyl.faces),1:3)=cell2mat([num2cell(cyl.faces,2)]);
            %             vtc(end,:)=z3ix;
            % %
            if vizz
                figure
                h=patch('vertices',cyl.vertices,'faces',cyl.faces,'facecolor','r','edgecolor','k');
            end
            [ncyl.vertices,ncyl.faces]=meshcheckrepair(cyl.vertices,cyl.faces,'meshfix');
            if vizz
                figure, patch('vertices',ncyl.vertices,'faces',ncyl.faces,'facecolor','r');
            end
            [ncyl.vertices,~,ncyl.faces]=s2m(ncyl.vertices,ncyl.faces,1,1,'tetgen');
            figure, patch('vertices',ncyl.vertices,'faces',ncyl.faces(:,1:3),'facecolor','r');

            ncyl.facevertexcdata=repmat(cyl.facevertexcdata(1,:),size(ncyl.vertices,1),1);
            %cyl=ncyl;
            %cyl.faces=cyl.faces(:,1:3);

            %             meshresample(h.Vertices,h.Faces,1);
            %           [nodecon,~,facecon]=s2m(cyl.vertices,[num2cell(cyl.faces,2);{flip(z3ix,2)}],0.2,1,'tetgen'); % generate a tetrahedral mesh of the cylinders

            %figure, patch(cyl,'facecolor','r')
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
            electrode.insulation(inscnt).vertices=nudgetomean(cyl.Vertices);
            electrode.insulation(inscnt).facevertexcdata=cyl.FaceVertexCData;
        catch
            electrode.insulation(inscnt).faces=cyl.faces;
            electrode.insulation(inscnt).vertices=nudgetomean(cyl.vertices);
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
            electrode.insulation(inscnt).vertices=nudgetomean(cyl.Vertices);
            electrode.insulation(inscnt).facevertexcdata=cyl.FaceVertexCData;
        catch
            electrode.insulation(inscnt).faces=cyl.faces;
            electrode.insulation(inscnt).vertices=nudgetomean(cyl.vertices);
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

save([ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname],'electrode');
% visualize
if vizz
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
        fv(cnt).vertices=electrode.insulation(ins).vertices;
        fv(cnt).faces=electrode.insulation(ins).faces;
        cnt=cnt+1;
    end
    for con=1:length(electrode.contacts)
        electrode.contacts(con).vertices=X*[electrode.contacts(con).vertices,ones(size(electrode.contacts(con).vertices,1),1)]';
        electrode.contacts(con).vertices=electrode.contacts(con).vertices(1:3,:)';
        elrender{side}(cnt)=patch(electrode.contacts(con));

        specsurf(elrender{side}(cnt),elspec.contact_color,aData);
        fv(cnt).vertices=electrode.contacts(con).vertices;
        fv(cnt).faces=electrode.contacts(con).faces;

        cnt=cnt+1;

    end

    axis equal
    view(0,0);
end





% export to .STL


fv=ea_concatfv(fv,0);
% fv=ea_mapcolvert2face(fv);
fv.faces=fv.faces(:,[3,2,1]);
ea_stlwrite(['bsc_vercise_directed.stl'],fv);

figure
plotmesh(fv.vertices,fv.faces)



% add meshel

% for con=1:length(electrode.contacts)
%     meshel.con{con}.faces=electrode.contacts(con).faces;
%     meshel.con{con}.vertices=round(electrode.contacts(con).vertices,50);
%     [meshel.con{con}.vertices,meshel.con{con}.faces]=meshcheckrepair(meshel.con{con}.vertices,meshel.con{con}.faces,'dup');
%     [meshel.con{con}.vertices,meshel.con{con}.faces]=meshcheckrepair(meshel.con{con}.vertices,meshel.con{con}.faces,'deep');
%
% end
%
% for ins=1:length(electrode.insulation)
%     meshel.ins{ins}.faces=electrode.insulation(ins).faces;
%     meshel.ins{ins}.vertices=round(electrode.insulation(ins).vertices,50);
%     [meshel.ins{ins}.vertices,meshel.ins{ins}.faces]=meshcheckrepair(meshel.ins{ins}.vertices,meshel.ins{ins}.faces,'dup');
%     [meshel.ins{ins}.vertices,meshel.ins{ins}.faces]=meshcheckrepair(meshel.ins{ins}.vertices,meshel.ins{ins}.faces,'deep');
% end






%ea_genvol_boston_dir(meshel,elspec,z3ix,vizz);



function node=nudgetomean(node)
nudge=0;


if nudge

    centroid=mean(node,1);
    tos=repmat(centroid,size(node,1),1)-node;
    tos=tos*0.01;
    node=node+tos;
end

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
