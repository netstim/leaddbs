function [emesh,nmesh,activeidx]=ea_mesh_electrode(fv,elfv,eltissuetype,electrode,options,S,side,elnumel)
% meshing an electrode and tissue structures bounded by a cylinder

%% load the nucleus data
tic
meshel=electrode.meshel;
vizz=0;
stlexport=1;
if vizz
    figure
    for f=1:length(fv)
        patch(fv(f),'FaceColor','none');
    end
    for f=1:length(elfv)
        patch(elfv(f),'FaceColor','none');
    end
    
%     figure
%     rgbcs='rgbcmyrgbcmyrgbcmy';
%     for f=1:length(elfv)
%         patch(elfv(f),'FaceColor',rgbcs(f),'FaceAlpha',0.5,'EdgeColor','none');
%     end
%     axis equal
%     
end

%% user defined parameters


orig=electrode.tail_position-3*(electrode.head_position-electrode.tail_position);
etop=electrode.head_position-3*(electrode.tail_position-electrode.head_position);
if vizz
    hold on
    plot3(orig(1),orig(2),orig(3),'y*');
    plot3(etop(1),etop(2),etop(3),'y*');
    oe=[orig;etop];
    plot3(oe(:,1),oe(:,2),oe(:,3),'r-');
end
%orig=[10.1,-15.07,-9.978];    % starting point of the electrode axis
%etop=[34.39,20.23,43.14];     % end point of the electrode axis
electrodelen=norm(etop-orig); % length of the electrode
electroderadius=0.65;            % radius of the electrode

electrodetrisize=0.2;  % the maximum triangle size of the electrode mesh
bcyltrisize=0.3;       % the maximum triangle size of the bounding cyl
nucleidecimate=0.2;    % downsample the nucleius mesh to 20%

cylz0=-35;     % define the lower end of the bounding cylinder
cylz1=35;     % define the upper end of the bounding cylinder
cylradius=15; % define the radius of the bounding cylinder  
ndiv=50;      % division of circle for the bounding cylinder

ncount=length(fv);     % the number of nuclei meshes inside fv()

v0=(etop-orig)/electrodelen;               % unitary dir
c0=[0 0 0];
v=[0 0 1];

%% loading the electrode surface model

ncyl=[];
fcyl=[];
scyl=[];
seeds=[];

for i=1:length(meshel.ins)
    fcyl=[fcyl; meshel.ins{i}.faces+size(ncyl,1)];
    if(i<length(meshel.ins))
        scyl=[scyl; meshel.ins{i}.endplates+size(ncyl,1)]; % had to rebuild the endplates
    end
    ncyl=[ncyl; meshel.ins{i}.vertices];
    seeds=[seeds; mean(meshel.ins{i}.vertices)];
end
for i=1:length(meshel.con)
    fcyl=[fcyl; meshel.con{i}.faces+size(ncyl,1)];
    scyl=[scyl; meshel.ins{i}.endplates+size(ncyl,1)];
    ncyl=[ncyl; meshel.con{i}.vertices];
    seeds=[seeds; mean(meshel.ins{i}.vertices)];
end

[ncyl, I, J]=unique(ncyl, 'rows');
fcyl=unique(round(J(fcyl)),'rows');
scyl=unique(round(J(scyl)),'rows');

fcyl=num2cell(fcyl,2);
scyl=num2cell(scyl,2);

%% convert the obtain the electrode surface mesh model
[node,elem,face]=s2m(ncyl,{fcyl{:}, scyl{:}},electrodetrisize,100,'tetgen',seeds,[]); % generate a tetrahedral mesh of the cylinders

%plotmesh(node,elem) % plot the electrode mesh for now

%% load the nucleus surfaces
nobj=[];
fobj=[];
nseeds=[];

for i=1:ncount
    no=fv(i).vertices;
    fo=fv(i).faces;
    [no,fo]=meshresample(no,fo,nucleidecimate); % mesh is too dense, reduce the density by 80%
    [no,fo]=meshcheckrepair(no,fo,'meshfix');  % clean topological defects
    fobj=[fobj;fo+size(nobj,1)];
    nobj=[nobj;no];
    nseeds=[nseeds; mean(no)];
end

%% merge the electrode mesh with the nucleus mesh
ISO2MESH_SURFBOOLEAN='cork';   % now intersect the electrode to the nucleus
[nboth,fboth]=surfboolean(node,face(:,[1 3 2]),'resolve',nobj,fobj);
clear ISO2MESH_SURFBOOLEAN;

%% create a bounding cylinder
%[anbcyl,afbcyl]=meshacylinder(orig, etop,cylradius,bcyltrisize,10,ndiv);


c0bbc=c0+cylz0*v;     
c1bbc=c0+cylz1*v;
[nbcyl,fbcyl]=meshacylinder(c0bbc, c1bbc,cylradius,bcyltrisize,10,ndiv);
nbcyl=rotatevec3d(nbcyl,v0,v);
nbcyl=nbcyl+repmat(orig,size(nbcyl,1),1);


%figure
%patch('vertices',anbcyl,'faces',afbcyl,'FaceColor','none','EdgeColor','b');
%patch('vertices',nbcyl,'faces',fbcyl,'FaceColor','none','EdgeColor','r');

seedbbc=nbcyl(1,:)*(1-1e-2)+mean(nbcyl)*1e-2;  % define a seed point for the bounding cylinder

%% cut the electrode+nucleus mesh by the bounding cylinder
ISO2MESH_SURFBOOLEAN='cork';
[nboth2,fboth2]=surfboolean(nbcyl,fbcyl(:,[1 3 2]),'first',nboth,fboth);
clear ISO2MESH_SURFBOOLEAN;

%% remove duplicated nodes in the surface
[nboth3,fboth3]=meshcheckrepair(nboth2,fboth2,'dup');
[nboth4,fboth4]=meshcheckrepair(nboth3,fboth3,'deep');

%% define seeds along the electrode axis
[t,baryu,baryv,faceidx]=raytrace(orig,v0,nboth4,fboth4);
t=sort(t(faceidx));
t=(t(1:end-1)+t(2:end))*0.5;
seedlen=length(t);
%electrodeseeds=repmat(orig(:)',seedlen,1)+repmat(v0(:)',seedlen,1).*repmat(t(:)-1,1,3);
%% create tetrahedral mesh of the final combined mesh (seeds are ignored, tetgen 1.5 automatically find regions)
[nmesh,emesh]=s2m(nboth3,fboth3,1,5);


%% remapping the region labels
etype=emesh(:,end);
labels=unique(etype);

eleccoord=rotatevec3d(nmesh-repmat(orig,size(nmesh,1),1),v,v0); % convert the centroids to the electrode cylinder coordinate

maxradius=zeros(length(labels),1);
zrange=zeros(length(labels),2);
centroids=zeros(length(labels),3);
for i=1:length(labels)
    centroids(i,:)=mean(meshcentroid(nmesh,emesh(etype==labels(i),1:3)));
    maxext(i,:)=max(nmesh(emesh(etype==labels(i),1),:));
    cc=(meshcentroid(eleccoord,emesh(etype==labels(i),1:3))); % centroids of each label
    maxradius(i)=sqrt(max(sum(cc(:,1:2).*cc(:,1:2),2)));
    zrange(i,:)=[min(cc(:,3)) max(cc(:,3))];
end


disp(['We have ',num2str(length(centroids)),' regions and need to map these to tissue types.']);
tissuelabels=zeros(length(centroids),1);
if vizz
%     h=figure;
%     plotmesh(nmesh,emesh(:,1:5),'linestyle','none','facealpha',0.1);
%     hold on
end


% init activeidx:

for s=1:4
    for c=1:elnumel
        activeidx(s).con(c).ix=[];
        activeidx(s).con(c).perc=0;
        activeidx(s).con(c).pol=0;
    end
end


active=find(S.activecontacts{side});

switch side
    case 1
        sidec='R';
    case 2
        sidec='L';
end
for reg=1:length(centroids)
   % first check if whether contact or insulator
   
   % a - check contacts:
    
   for con=find(eltissuetype==3);
        convin=ea_intriangulation(elfv(con).vertices,elfv(con).faces,centroids(reg,:));
        thiscompsnodes=emesh(emesh(1:end,5)==reg,1:4); % get this components nodes
        dirinodes=nmesh(thiscompsnodes,:);
        dirinodes=ea_nudgedirinodes(dirinodes,centroids(reg,:));
        in=double(ea_intriangulation(elfv(con).vertices,elfv(con).faces,dirinodes));

        if convin && mean(in)>0.7
            if ismember(con,active)
                
                % we captured an active contact. need to assign to correct
                % source and polarity
                for source=1:4
                   if S.([sidec,'s',num2str(source)]).amp % then this active contact could be from this source since source is active
                       if S.([sidec,'s',num2str(source)]).(['k',num2str(con-1)]).perc % current captured contact is from this source    
                           activeidx(source).con(con).ix=[activeidx(source).con(con).ix;unique(thiscompsnodes(:))];
                           activeidx(source).con(con).pol=S.([sidec,'s',num2str(source)]).(['k',num2str(con-1)]).pol;
                           activeidx(source).con(con).perc=S.([sidec,'s',num2str(source)]).(['k',num2str(con-1)]).perc;
                       end
                   end
                end
                
                
            end
            tissuelabels(reg)=3; % set contact
            disp(['Region ',num2str(reg),' captured by contact material.']);
            if vizz
%                  figure('name',['Conducting region ',num2str(reg)]);
%                  hold on
%                  patch('vertices',elfv(con).vertices,'faces',elfv(con).faces,'FaceColor','none','EdgeColor','b');
%                  patch('vertices',nmesh,'faces',emesh(emesh(:,5)==reg,1:4),'FaceColor','none','EdgeColor','r');
%                  plot3(centroids(reg,1),centroids(reg,2),centroids(reg,3),'go');
%                  axis equal
            end
            break
        end
    end
   if tissuelabels(reg); continue; end % move to next component if already assigned.
 
   % b - check insulation:
    
   for ins=find(eltissuetype==4);
        convin=ea_intriangulation(elfv(ins).vertices,elfv(ins).faces,centroids(reg,:));
        dirinodes=nmesh(emesh(emesh(1:end,5)==reg,1:4),:);
        dirinodes=ea_nudgedirinodes(dirinodes,centroids(reg,:));
        in=double(ea_intriangulation(elfv(ins).vertices,elfv(ins).faces,dirinodes));
        if convin && mean(in)>0.7
            tissuelabels(reg)=4; % set insulation
            disp(['Region ',num2str(reg),' captured by insulating material.']);
            if vizz
%                 figure('name',['Insulating region ',num2str(reg)]);
%                 hold on
%                 patch('vertices',elfv(ins).vertices,'faces',elfv(ins).faces,'FaceColor','none','EdgeColor','b');
%                 patch('vertices',nmesh,'faces',emesh(emesh(:,5)==reg,1:4),'FaceColor','none','EdgeColor','r');
%                 plot3(centroids(reg,1),centroids(reg,2),centroids(reg,3),'go');
%                 axis equal
            end
            break
        end
    end
   if tissuelabels(reg); continue; end % move to next component if already assigned.
     
    
    % if not: if grey matter, then white matter
    
    for gm=1:length(fv)
        convin=ea_intriangulation(fv(gm).vertices,fv(gm).faces,centroids(reg,:));
        dirinodes=nmesh(emesh(emesh(1:end,5)==reg,1:4),:);
        dirinodes=ea_nudgedirinodes(dirinodes,centroids(reg,:));
        in=double(ea_intriangulation(fv(gm).vertices,fv(gm).faces,dirinodes));
       
        if convin && mean(in)>0.7
            tissuelabels(reg)=1; % set grey matter
            disp(['Region ',num2str(reg),' captured by grey matter.']);
            if vizz
%                 figure('name',['GM Region ',num2str(reg)]);
%                 hold on
%                 patch('vertices',fv(gm).vertices,'faces',fv(gm).faces,'FaceColor','none','EdgeColor','b');
%                 patch('vertices',nmesh,'faces',emesh(emesh(:,5)==reg,1:4),'FaceColor','none','EdgeColor','r');
%                 plot3(centroids(reg,1),centroids(reg,2),centroids(reg,3),'go');
%                 axis equal
            end
            break
        end
    end
    if tissuelabels(reg); continue; end
    
    
    % assign the rest to white matter
    tissuelabels(reg)=2; % set white matter
    
    if vizz
%         figure('name',['WM Region ',num2str(reg)]);
%         hold on
%         patch('vertices',nmesh,'faces',emesh(emesh(:,5)==reg,1:4),'FaceColor','none','EdgeColor','r');
%         for g=1:length(fv)
%         patch('vertices',fv(g).vertices,'faces',fv(g).faces,'FaceColor','none','EdgeColor','b');
%         end
%         for e=1:length(elfv)
%         patch('vertices',elfv(e).vertices,'faces',elfv(e).faces,'FaceColor','none','EdgeColor','b'); 
%         end
%         plot3(centroids(reg,1),centroids(reg,2),centroids(reg,3),'go');
%         axis equal
    end
    
end


%gmlabels=setdiff(labels,[wmlabels; electrodelabel]); % the remaining ones are from nuclei meshes.

tissuetype=emesh(:,5);
for tt=1:4
    tl=find(tissuelabels==tt);
    tissuetype(ismember(emesh(:,5),tl))=tt;
    disp([num2str(length(tl)),' components with ',num2str(sum(tissuetype==tt)),' tetraeders total assigned to tissue type ',num2str(tt),'.']);
end

emesh(:,5)=tissuetype;

if vizz
%% plot the final tetrahedral mesh
figure
hold on;
plotmesh(nmesh,emesh,'linestyle','none','facealpha',0.2)
end


if stlexport
    tissuelabels={'grey','white','contacts','insulation'};
    mkdir([options.root,options.patientname,filesep,'headmodel',filesep]);
    for tt=1:length(tissuelabels)
        savestl(nmesh,emesh(emesh(:,5)==tt,1:4),[options.root,options.patientname,filesep,'headmodel',filesep,tissuelabels{tt},'.stl'],tissuelabels{tt});
    end
end
% plot all 4 tissue types:
if vizz
    for t=1:4
        figure('Name',['Tissue ',num2str(t)]);
        hold on;
        plotmesh(nmesh,emesh(emesh(:,5)==t,:),'linestyle','none','facealpha',0.2)
    end
end


toc



function nodes=ea_nudgedirinodes(nodes,centroid)
% get to ~1000 comparison points
div=round(length(nodes)/1000);
if div<1
    div=1;
end
nodes=nodes(1:div:end,:);
nodes=nodes.*9;
nodes=nodes+repmat(centroid,length(nodes),1);
nodes=nodes./10;

