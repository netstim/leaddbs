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
electrodetrisize=0.1;
options.sides=1;
elstruct.name=options.elmodel;
options=ea_resolve_elspec(options);
elspec=options.elspec;
resultfig=figure('visible','off');

ncnt=1; % count for numeshel
concnt=1;
inscnt=1;
jetlist=jet;
%   jetlist=jet;
N=60; % resolution of electrode points
aData=1;

cnt=4;

scf=1; % scale factor for precision problems, do not modify.
elspec.contact_length=elspec.contact_length*scf;
elspec.contact_diameter=elspec.contact_diameter*scf;
elspec.contact_spacing=elspec.contact_spacing*scf;
elspec.lead_diameter=elspec.lead_diameter*scf;
elspec.tip_diameter=elspec.tip_diameter*scf;
elspec.tip_length=elspec.tip_length*scf;

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
    %[cX,cY,cZ] = cylinder((diams),(2*N-1));
    
    
    %figure, plot3(cX,cY,cZ,'r.')
    
    cZ=cZ.*(startpoint(3)-lowerpoint(3)); % scale to fit tip-diameter
    cZ=cZ+lowerpoint(3);
    
    p=surf2patch(surf(cX,cY,cZ),'triangles');
    
    
    % add meshing-version to shaft
    %[cX,cY,cZ] = ea_singlecylinder((diams),N);
    [cX,cY,cZ] = cylinder((diams),(2*N-1));
    
    cZ=cZ.*(startpoint(3)-lowerpoint(3)); % scale to fit tip-diameter
    cZ=cZ+lowerpoint(3);
    a=surf2patch(surf(cX,cY,cZ));
    a=ea_reordercylinder(a,2);
    
    meshel.ins{1}.faces=a.faces;
    meshel.ins{1}.vertices=a.vertices;
    ndiv=length(meshel.ins{1}.vertices)/2;
    meshel.ins{1}.endplates=[1:ndiv;ndiv+1:2*ndiv];
    
    % nbefore we had two endplates for each, now only top:
    %    numeshel(ncnt).faces=[num2cell(a.faces,2);num2cell([1:ndiv;ndiv+1:2*ndiv],2)];
    
    
    numeshel(ncnt).faces=[num2cell(a.faces,2);num2cell([1:ndiv;ndiv+1:2*ndiv],2)];
    numeshel(ncnt).vertices=[a.vertices];
    eltype(ncnt)=2; % ins
    
    [electrode.insulation(inscnt).vertices,~,electrode.insulation(inscnt).faces]=s2m(numeshel(ncnt).vertices,{numeshel(ncnt).faces{:}},electrodetrisize,100,'tetgen',[],[]); % generate a tetrahedral mesh of the cylinders
    electrode.insulation(inscnt).faces=electrode.insulation(inscnt).faces(:,1:3);
    inscnt=inscnt+1;
    

    
    ncnt=ncnt+1;
    
    
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
    for cntct=1:4 % first contact is the tip (see below).
        
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
            
            %[cX,cY,cZ] = ea_singlecylinder((tipdiams),N);
            [cX,cY,cZ] = cylinder((tipdiams),(2*N-1));
            
            cZ=cZ.*(elspec.tip_length); % scale to fit tip-diameter
            
            
            
            cX=cX+X1(1);
            cY=cY+X1(2);
            cZ=cZ-(2*elspec.tip_length)/2+X1(3);
            a=surf2patch(surf(cX,cY,cZ));
            
            a=ea_reordercylinder(a,2);
            meshel.con{1}.faces=a.faces;
            meshel.con{1}.vertices=a.vertices;
            ndiv=length(meshel.con{1}.vertices)/12;
            meshel.con{1}.endplates=[find(a.vertices(:,3)==max(a.vertices(:,3)))'];
            
            %             figure
            %             patch(a,'FaceColor','none');
            %             hold on
            %             patch('vertices',a.vertices,'faces',find(a.vertices(:,3)==max(a.vertices(:,3)))','facecolor','r')
            %             keyboard
            numeshel(ncnt).faces=[num2cell(a.faces,2);num2cell([find(a.vertices(:,3)==max(a.vertices(:,3)))'],2)];
            numeshel(ncnt).vertices=a.vertices;
            eltype(ncnt)=1; % con
            
            [electrode.contacts(concnt).vertices,~,electrode.contacts(concnt).faces]=s2m(numeshel(ncnt).vertices,{numeshel(ncnt).faces{:}},electrodetrisize,100,'tetgen',[],[]); % generate a tetrahedral mesh of the cylinders
            electrode.contacts(concnt).faces=electrode.contacts(concnt).faces(:,1:3);
            concnt=concnt+1;
            
            
            %             figure
            %             hold on
            %             for c=1:length(numeshel(ncnt).faces)
            %             patch('vertices',numeshel(ncnt).vertices,'faces',numeshel(ncnt).faces{c},'facecolor',rand(3,1),'facealpha',0.3);
            %             end
            %
            ncnt=ncnt+1;
            
            
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
            
            
        elseif cntct==2 || cntct==3 % segmented cylinders
            % define two points to define cylinder.
            X1=coords_mm{side}(cntct,:)+trajvector*(elspec.contact_length/2);
            [no,fc,seeds,fcc,noc,fci,noi] = ea_segmented_cylinder_qq_ah(60,1,0.75,1,3,0.8);
            
            %     figure, plotmesh(no,fcc);
            %     figure, plotmesh(no,fci);
    
            %[ino,fci]=removedupnodes(ino,fci,1e-6);
            %[cno,fci]=removedupnodes(cno,fcc,1e-6);
%            [cno,fci]=removeisolatednode(no,fcc);

 noc(:,1)=noc(:,1).*(elspec.contact_diameter/2);
            noc(:,2)=noc(:,2).*(elspec.contact_diameter/2);
            noc(:,3)=noc(:,3).*(elspec.contact_length);
            
            noc(:,1)=noc(:,1)+X1(1);
            noc(:,2)=noc(:,2)+X1(2);
            noc(:,3)=noc(:,3)+X1(3);


 noi(:,1)=noi(:,1).*(elspec.contact_diameter/2);
            noi(:,2)=noi(:,2).*(elspec.contact_diameter/2);
            noi(:,3)=noi(:,3).*(elspec.contact_length);
            
            noi(:,1)=noi(:,1)+X1(1);
            noi(:,2)=noi(:,2)+X1(2);
            noi(:,3)=noi(:,3)+X1(3);

electrode.contacts(concnt).vertices=noc;
electrode.contacts(concnt).faces=fcc(:,1:3);
            concnt=concnt+1;
            
electrode.insulation(inscnt).vertices=noi;
electrode.insulation(inscnt).faces=fci(:,1:3);
            inscnt=inscnt+1;
            
                
              
            
            % shape up:
            [vnode,velem,vface]=s2m(no,fc,1,3);
            
            
            % %% VARIANT B: Use cell format and use non-voled version first
            no(:,1)=no(:,1).*(elspec.contact_diameter/2);
            no(:,2)=no(:,2).*(elspec.contact_diameter/2);
            no(:,3)=no(:,3).*(elspec.contact_length);
            
            no(:,1)=no(:,1)+X1(1);
            no(:,2)=no(:,2)+X1(2);
            no(:,3)=no(:,3)+X1(3);
            
            [volnode,volelem,volface]=s2m(no,fc,1,3);
%             figure
%             subplot(121);
%             plotmesh(volnode,volelem);
            



            % Here we feed the whole segmented cylinder - including
            % insulation - to the .con element of meshel since division is
            % not needed at this point.
            
            numeshel(ncnt).faces=fc;
            numeshel(ncnt).vertices=[no];
            eltype(ncnt)=3; % mixed.
            %             figure
            %             hold on
            %             for c=1:length(numeshel(ncnt).faces)
            %             patch('vertices',numeshel(ncnt).vertices,'faces',numeshel(ncnt).faces{c},'facecolor',rand(3,1),'facealpha',0.3);
            %             end
            %
            
            ncnt=ncnt+1;
            %
            %     figure,
            %     plot3(numeshel(ncnt).vertices(:,1),numeshel(ncnt).vertices(:,2),numeshel(ncnt).vertices(:,3),'r.')
            %
            %     hold on
            %     ncnt=3;
            %     plot3(numeshel(ncnt).vertices(:,1),numeshel(ncnt).vertices(:,2),numeshel(ncnt).vertices(:,3),'r.')
            
            vcnt=1;
            for contact=1:3 % ring segment
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
                
                fvseg(contact).faces=elrender{side}(cnt).Faces;
                fvseg(contact).vertices=elrender{side}(cnt).Vertices;
                
                if vizz
                    %figure
                    %hold on
                    %plotmesh(meshel.con{end}.vertices,meshel.con{end}.faces)
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
                patch('Vertices',thisno,'Faces',thisfc,'FaceColor','none','EdgeColor',rand(3,1));
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
                
                % meshing version, these contacts are added to meshel.con
                % including their insulating space since doesn't matter.
                
                meshel.ins{end+1}.faces=elrender{side}(cnt).Faces;
                meshel.ins{end}.vertices=elrender{side}(cnt).Vertices;
                
                fvseg(ins).faces=elrender{side}(cnt).Faces;
                fvseg(ins).vertices=elrender{side}(cnt).Vertices;
                %figure, patch('vertices',fvseg(ins).vertices,'faces',fvseg(ins).faces,'facecolor','none');

                
                log(cnt)=0; % insulation
                
                cnt=cnt+1;
                vcnt=vcnt+1;
                
            end
            
            
            
            
            %             plotmesh(contno,contfc);
            %             subplot(122);
            %             plotmesh(insuno,insufc);
            %
            %             figure
            %             plotmesh(vnode,vface)
            
            
            
            
            
            
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
            %[cX,cY,cZ] = ea_singlecylinder((diams),N);
            [cX,cY,cZ] = cylinder((diams),(2*N-1));
            
            cZ=cZ.*(elspec.contact_length); % scale to fit tip-diameter
            htd=(max(cZ(:))/2);
            cZ=cZ-htd;
            cZ=cZ+coords_mm{side}(cntct,3);
            a=surf2patch(surf(cX,cY,cZ));
            
            
            
            
            a=ea_reordercylinder(a,2);
            meshel.con{end+1}.faces=a.faces;
            meshel.con{end}.vertices=a.vertices;
            ndiv=length(meshel.con{end}.vertices)/2;
            %             figure
            %             patch(a,'FaceColor','none')
            %             hold on
            %             patch('vertices',a.vertices,'faces',1:ndiv,'facecolor','r');
            %             patch('vertices',a.vertices,'faces',ndiv+1:2*ndiv,'facecolor','b');
            
            meshel.con{end}.endplates=[1:ndiv;ndiv+1:2*ndiv];
            
            
            
            numeshel(ncnt).faces=[num2cell(a.faces,2);num2cell([1:ndiv;ndiv+1:2*ndiv],2)];
            numeshel(ncnt).vertices=a.vertices;
            eltype(ncnt)=1; % con
            %             figure
            %             hold on
            %             for c=1:length(numeshel(ncnt).faces)
            %             patch('vertices',numeshel(ncnt).vertices,'faces',numeshel(ncnt).faces{c},'facecolor',rand(3,1));
            %             end
            
            [electrode.contacts(concnt).vertices,~,electrode.contacts(concnt).faces]=s2m(numeshel(ncnt).vertices,{numeshel(ncnt).faces{:}},electrodetrisize,100,'tetgen',[],[]); % generate a tetrahedral mesh of the cylinders
            electrode.contacts(concnt).faces=electrode.contacts(concnt).faces(:,1:3);
            concnt=concnt+1;
            
            ncnt=ncnt+1;
            
            
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
        %[cX,cY,cZ] = ea_singlecylinder((diams),N);
        [cX,cY,cZ] = cylinder((diams),(2*N-1));
        
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
        
        %if cntct==3
        numeshel(ncnt).faces=[num2cell(a.faces,2);num2cell([1:ndiv;ndiv+1:2*ndiv],2)]; % only need an top endplate in case of last isolating segment
        %else
        %    numeshel(ncnt).faces=[num2cell(a.faces,2)];
        %end
        numeshel(ncnt).vertices=a.vertices;
        eltype(ncnt)=2; % ins
        
        [electrode.insulation(inscnt).vertices,~,electrode.insulation(inscnt).faces]=s2m(numeshel(ncnt).vertices,{numeshel(ncnt).faces{:}},electrodetrisize,100,'tetgen',[],[]); % generate a tetrahedral mesh of the cylinders
        electrode.insulation(inscnt).faces=electrode.insulation(inscnt).faces(:,1:3);
        inscnt=inscnt+1;
        
        %         figure
        %             hold on
        %             for c=1:length(numeshel(ncnt).faces)
        %             patch('vertices',numeshel(ncnt).vertices,'faces',numeshel(ncnt).faces{c},'facecolor',rand(3,1));
        %             end
        %         keyboard
        ncnt=ncnt+1;
        
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



% export vol version:

ea_genvol_bscdirected(numeshel,elspec,1);


% more stuff for the visualization version (electrode.contacts &
% electrodes.insulation)


electrode.electrode_model=elstruct.name;
electrode.head_position=[0,0,0.5*options.elspec.contact_length];
electrode.tail_position=[0,0,4*options.elspec.contact_length+(4-1)*options.elspec.contact_spacing-0.5*options.elspec.contact_length];

electrode.x_position=[options.elspec.lead_diameter/2,0,0.5*options.elspec.contact_length];
electrode.y_position=[0,options.elspec.lead_diameter/2,0.5*options.elspec.contact_length];

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

%figure
%plotmesh(fv.vertices,fv.faces)



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
