function [oemesh,nmesh,activeidx,wmboundary,centroids,tissuetype]=ea_mesh_electrode(fv,elfv,eltissuetype,electrode,options,S,side,elnumel,transformmatrix,elspec)
% meshing an electrode and tissue structures bounded by a cylinder
% load the nucleus data
ea_dispt('Generating tetraedrical mesh...');
% meshel=electrode.meshel;
vizz=0;
stlexport=0;
if vizz
    figure
    for f=1:length(fv)
        patch(fv(f),'FaceColor','none');
    end
    for f=1:length(elfv)
        patch(elfv(f),'FaceColor','none');
    end
end

if max(S.amplitude{side})>4
    stretchfactor= 0.3125*(max(S.amplitude{side})/2.5);
else
    stretchfactor=0.5;
end
%stretchfactor=0.5*(max(S.amplitude{side})/2.5);

orig=electrode.tail_position-3*stretchfactor*(electrode.head_position-electrode.tail_position);
etop=electrode.head_position-3*stretchfactor*(electrode.tail_position-electrode.head_position);

el_o_orig=[0,0,15+(20*stretchfactor)];
el_o_etop=[0,0,-20*stretchfactor];

nucleidecimate=0.2;    % downsample the nucleius mesh to 20%

bcyltrisize=0.01;       % the maximum triangle size of the bounding cyl

cylz0=-30;     % define the upper end of the bounding cylinder
cylz1=30;     % define the lowe end of the bounding cylinder
cylradius=40*stretchfactor; % define the radius of the bounding cylinder

ndiv=50;      % division of circle for the bounding cylinder
electrodelen=norm(etop-orig); % length of the electrode

v0=(etop-orig)/electrodelen;               % unitary dir
c0=[0 0 0];
v=[0 0 1];

elmodel_fn=[ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname,'_vol.mat'];
if ~exist(elmodel_fn,'file')
    ea_generate_electrode_specs; % regenerate all electrode specifications
end

elmodel = load(elmodel_fn);
node = elmodel.node;
face = elmodel.face;

% apply transformation matrix to electrode nodes:
node=transformmatrix*[node,ones(size(node,1),1)]';
node=node(1:3,:)';

% - this is the node / elem / face made by tetgen of the electrode only.
if vizz
    figure
    fvv.faces=face(:,1:3);
    fvv.vertices=node;
    patch(fvv,'edgecolor','b','facecolor','none');
    axis equal
end

%plotmesh(node,elem) % plot the electrode mesh for now

% create a bounding cylinder
%[anbcyl,afbcyl]=meshacylinder(orig, etop,cylradius,bcyltrisize,10,ndiv);

c0bbc=c0+cylz0*v;
c1bbc=c0+cylz1*v;
%     [nbcyl,fbcyl]=meshacylinder(c0bbc, c1bbc,cylradius,bcyltrisize,10,ndiv);
%     nbcyl=rotatevec3d(nbcyl,v0,v);
%     nbcyl=nbcyl+(repmat(orig,size(nbcyl,1),1)/stretchfactor);

[nbcyl,fbcyl]=meshacylinder(el_o_etop,el_o_orig,cylradius,bcyltrisize,10,ndiv);
nbcyl=transformmatrix*[nbcyl,ones(length(nbcyl),1)]';
nbcyl=nbcyl(1:3,:)';

if vizz
    figure
    fva.faces=fbcyl(:,1:3);
    fva.vertices=nbcyl;
    patch(fva,'edgecolor','m','facecolor','none');
    hold on
    plot3(orig(1),orig(2),orig(3),'r*');
    plot3(etop(1),etop(2),etop(3),'g*');
    plot3(electrode.tail_position(1),electrode.tail_position(2),electrode.tail_position(3),'k*');
    plot3(electrode.head_position(1),electrode.head_position(2),electrode.head_position(3),'b*');
    axis equal
end

% if isempty(fv) % use TPM
%     c1=ea_load_nii([ea_space(options),'TPM.nii,1']);
%     voxnbcyl=c1.mat\[nbcyl,ones(length(nbcyl),1)]';
%     voxnbcyl=voxnbcyl(1:3,:)';
%     cyl=surf2vol(voxnbcyl,fbcyl,1:size(c1.img,2),1:size(c1.img,1),1:size(c1.img,3));
%     cyl=imfill(cyl,'holes');
%
%     cyl=double(smooth3(cyl,'gaussian',[3 3 3]));
%     c1.img=c1.img.*permute(cyl,[2,1,3]);
%     fv=isosurface(c1.img,0.5,'noshare');
%     fv.vertices=c1.mat*[fv.vertices,ones(length(fv.vertices),1)]';
%     fv.vertices=fv.vertices(1:3,:)';
%     tpmuse=1;
% else
tpmuse=0;
% end

if ~isempty(fv) % use atlas to define GM
    % load the nucleus surfaces
    nobj=[];
    fobj=[];
    ncount=length(fv);     % the number of nuclei meshes inside fv()
    ISO2MESH_SURFBOOLEAN='cork';   % now intersect the electrode to the nucleus

    if ncount==1 && isempty(fv(1).vertices) % no gray matter
        graymatterpresent=0;
    else
        graymatterpresent=1;
        for i=1:ncount
            no=fv(i).vertices;
            fo=fv(i).faces;
            [no,fo]=meshresample(no,fo,nucleidecimate); % mesh is too dense, reduce the density by 80%
            % [no,fo]=meshcheckrepair(no,fo,'meshfix');  % clean topological defects

            % merge all nuclei
            if isempty(nobj)
                nobj=no;
                fobj=fo;
            else
                [nobj,fobj]=surfboolean(no,fo,'resolve',nobj,fobj);
            end
            % fobj=[fobj;fo+size(nobj,1)];
            % nobj=[nobj;no];
        end

        if vizz
            figure
            patch('Vertices',nobj,'Faces',fobj,'FaceColor','none');
            patch('Vertices',node,'Faces',face(:,1:3),'FaceColor','blue');
        end
    end
    % merge the electrode mesh with the nucleus mesh
else
    graymatterpresent=0;
end

if graymatterpresent
    [nboth,fboth]=surfboolean(node,face(:,1:3),'resolve',nobj,fobj);
else
    nboth=node;
    fboth=face;
end

clear ISO2MESH_SURFBOOLEAN;

if vizz
    figure
    fvv.faces=fboth(:,1:3);
    fvv.vertices=nboth;
    patch(fvv,'edgecolor','b','facecolor','none');
end

%figure
%patch('vertices',anbcyl,'faces',afbcyl,'FaceColor','none','EdgeColor','b');
%patch('vertices',nbcyl,'faces',fbcyl,'FaceColor','none','EdgeColor','r');

%seedbbc=nbcyl(1,:)*(1-1e-2)+mean(nbcyl)*1e-2;  % define a seed point for the bounding cylinder

% cut the electrode+nucleus mesh by the bounding cylinder
ISO2MESH_SURFBOOLEAN='cork';

[nboth2,fboth2]=surfboolean(nbcyl,fbcyl(:,[1 3 2]),'resolve',nboth,fboth);
%[nboth2,fboth2]=surfboolean(nbcyl,fbcyl(:,[1 3 2]),'first',nboth2,fboth2);

clear ISO2MESH_SURFBOOLEAN;
if vizz
    figure('name','nboth2');
    %patch('vertices',nbothbc,'faces',fbothbc,'FaceColor','none')
    fvv.faces=fboth2(:,1:3);
    fvv.vertices=nboth2;
    patch(fvv,'edgecolor','green','facecolor','none');
    axis equal
end

% remove duplicated nodes in the surface
[nboth3,fboth3]=meshcheckrepair(nboth2,fboth2,'dup');
[nboth3,fboth3]=meshcheckrepair(nboth3,fboth3,'deep');

%figure, patch('faces',fboth4,'vertices',nboth4,'facecolor','r','facealpha',0.3);
if vizz
    figure('name','nboth3');
    fvv.faces=fboth3(:,1:3);
    fvv.vertices=nboth3;
    patch(fvv,'edgecolor','m','facecolor','none');
    axis equal
end

% define seeds along the electrode axis
%[t,baryu,baryv,faceidx]=raytrace(orig,v0,nboth4,fboth4);
%t=sort(t(faceidx));
%t=(t(1:end-1)+t(2:end))*0.5;
%seedlen=length(t);
%electrodeseeds=repmat(orig(:)',seedlen,1)+repmat(v0(:)',seedlen,1).*repmat(t(:)-1,1,3);

% create tetrahedral mesh of the final combined mesh (seeds are ignored, tetgen 1.5 automatically find regions)
% - this is the part where we have all 4 element types combined already.

[nmesh,emesh,fmesh]=s2m(nboth3,fboth3,1,3);
if vizz
    figure('name','Final mesh');
    fvv.faces=face(:,1:3);
    fvv.vertices=nmesh;
    patch(fvv,'edgecolor','k','facecolor','none');
    axis equal
end

% remapping the region labels
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
    % h=figure;
    % plotmesh(nmesh,emesh(:,1:5),'linestyle','none','facealpha',0.1);
    % hold on
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

wmboundary=[];

if vizz
    h=figure;
end

for reg=1:length(centroids)
    % first check if whether contact or insulator
    thiscompsnodes=emesh(emesh(1:end,5)==reg,1:4); % get this components nodes
    Ntc=size(thiscompsnodes,1);
    if Ntc>2500 % take a representative sample if whole points too large.
        thiscompsnodes=thiscompsnodes(round(linspace(1,Ntc,2500)),:);
    end
    tetrcents=mean(cat(3,nmesh(thiscompsnodes(:,1),:),nmesh(thiscompsnodes(:,2),:),nmesh(thiscompsnodes(:,3),:),nmesh(thiscompsnodes(:,4),:)),3);

    % a - check contacts:
    for con=find(eltissuetype==3)
        in=double(ea_intriangulation(elfv(con).vertices,elfv(con).faces,tetrcents));

        if vizz
            set(h,'name',num2str(mean(in)));
            hold off
            plot3(elfv(con).vertices(:,1),elfv(con).vertices(:,2),elfv(con).vertices(:,3),'r*');
            hold on
            plot3(tetrcents(:,1),tetrcents(:,2),tetrcents(:,3),'g.');
            drawnow
        end

        if (mean(in)>0.7)
            if ismember(con,active)
                % we captured an active contact. need to assign to correct
                % source and polarity
                for source=1:4
                    if S.([sidec,'s',num2str(source)]).amp % then this active contact could be from this source since source is active
                        if S.([sidec,'s',num2str(source)]).(['k',num2str(con+8*(side-1)-1)]).perc % current captured contact is from this source
                            activeidx(source).con(con).ix=[activeidx(source).con(con).ix;unique(thiscompsnodes(:))];
                            activeidx(source).con(con).pol=S.([sidec,'s',num2str(source)]).(['k',num2str(con+8*(side-1)-1)]).pol;
                            activeidx(source).con(con).perc=S.([sidec,'s',num2str(source)]).(['k',num2str(con+8*(side-1)-1)]).perc;
                        end
                    end
                end
            end
            tissuelabels(reg)=3; % set contact
            disp(['Region ',num2str(reg),' captured by contact material.']);
            if vizz
                % figure('name',['Conducting region ',num2str(reg)]);
                % hold on
                % patch('vertices',elfv(con).vertices,'faces',elfv(con).faces,'FaceColor','none','EdgeColor','b');
                % patch('vertices',nmesh,'faces',emesh(emesh(:,5)==reg,1:4),'FaceColor','none','EdgeColor','r');
                % plot3(centroids(reg,1),centroids(reg,2),centroids(reg,3),'go');
                % axis equal
            end
            break
        end
    end

    if tissuelabels(reg); continue; end % move to next component if already assigned.

    % b - check insulation:
    for ins=find(eltissuetype==4)
        in=double(ea_intriangulation(elfv(ins).vertices,elfv(ins).faces,tetrcents));
        if vizz
            set(h,'name',num2str(mean(in)));
            hold off
            plot3(elfv(ins).vertices(:,1),elfv(ins).vertices(:,2),elfv(ins).vertices(:,3),'r*');
            hold on
            plot3(tetrcents(:,1),tetrcents(:,2),tetrcents(:,3),'g.');
            drawnow
        end

        if (mean(in)>0.7)
            tissuelabels(reg)=4; % set insulation
            disp(['Region ',num2str(reg),' captured by insulating material.']);
            if vizz
                % figure('name',['Insulating region ',num2str(reg)]);
                % hold on
                % patch('vertices',elfv(ins).vertices,'faces',elfv(ins).faces,'FaceColor','none','EdgeColor','b');
                % patch('vertices',nmesh,'faces',emesh(emesh(:,5)==reg,1:4),'FaceColor','none','EdgeColor','r');
                % plot3(centroids(reg,1),centroids(reg,2),centroids(reg,3),'go');
                % axis equal
            end
            break
        end
    end

    if tissuelabels(reg); continue; end % move to next component if already assigned.

    % if not: if grey matter, then white matter
    if ~tpmuse
        if graymatterpresent
            for gm=1:length(fv)

                in=double(ea_intriangulation(fv(gm).vertices,fv(gm).faces,tetrcents));
                if vizz
                    set(h,'name',num2str(mean(in)));
                    hold off
                    plot3(fv(gm).vertices(:,1),fv(gm).vertices(:,2),fv(gm).vertices(:,3),'r*');
                    hold on
                    plot3(tetrcents(:,1),tetrcents(:,2),tetrcents(:,3),'g.');
                    drawnow
                end

                if (mean(in)>0.7)
                    tissuelabels(reg)=1; % set grey matter
                    disp(['Region ',num2str(reg),' captured by grey matter.']);
                    if vizz
                        figure('name',['GM Region ',num2str(reg)]);
                        hold on
                        patch('vertices',fv(gm).vertices,'faces',fv(gm).faces,'FaceColor','none','EdgeColor','b');
                        patch('vertices',nmesh,'faces',emesh(emesh(:,5)==reg,1:4),'FaceColor','none','EdgeColor','r');
                        plot3(centroids(reg,1),centroids(reg,2),centroids(reg,3),'go');
                        axis equal
                    end
                    break
                end
            end
        end
    else
        thisc=c1.mat\[centroids(reg,:),1]';
        pval=spm_sample_vol(c1,thisc(1),thisc(2),thisc(3),1);
        if pval>0.5 % GM
            tissuelabels(reg)=1;
            disp(['Region ',num2str(reg),' captured by grey matter.']);
        end
    end

    if tissuelabels(reg); continue; end

    % assign the rest to white matter: (this following code will not be executed if
    % label has already been assigned above).
    tissuelabels(reg)=2; % set white matter
    disp(['Region ',num2str(reg),' captured by white matter.']);
end

% now we need to get surface nodes based on nbcyl:
if isempty(which('rangesearch'))
    ea_error('Matlab Statistics Toolbox not installed. This is (unfortunately) needed to calculate VTAs this way.');
end
wmboundary=rangesearch(nmesh,nbcyl,0.1);

wmboundary=unique(cell2mat(wmboundary'));

if vizz
    figure,
    hold on
    plot3(nmesh(wmboundary,1),nmesh(wmboundary,2),nmesh(wmboundary,3),'r.');
    plot3(nbcyl(:,1),nbcyl(:,2),nbcyl(:,3),'b.');
    plot3(nmesh(:,1),nmesh(:,2),nmesh(:,3),'g.');
end

% [nbothbc]=surfboolean(nmesh,fmesh(:,1:3),'second',nbcyl,fbcyl);
% [~,wmboundary]=ismember(nbothbc,nmesh,'rows');
% wmboundary(wmboundary==0)=[];
% gmlabels=setdiff(labels,[wmlabels; electrodelabel]); % the remaining ones are from nuclei meshes.

tissuetype=emesh(:,5);
for tt=1:4
    tl=find(tissuelabels==tt);
    tissuetype(ismember(emesh(:,5),tl))=tt;
    disp([num2str(length(tl)),' components with ',num2str(sum(tissuetype==tt)),' tetraeders total assigned to tissue type ',num2str(tt),'.']);
end

if ~all(tissuetype)
    tissuetype(tissuetype==0)=2; % assign misfit tetraeders to white matter.
end
oemesh=emesh;
emesh(:,5)=tissuetype;

if vizz
    % plot the final tetrahedral mesh
    figure
    hold on;
    plotmesh(nmesh,emesh,'linestyle','none','facealpha',0.2)
end

if stlexport
    ea_dispt('Exporting STL files');
    tissuelabels={'grey','white','contacts','insulation'};

    headmodelDir = fullfile(options.subj.subjDir, 'headmodel', ea_nt(options));
    ea_mkdir(headmodelDir);
    filePrefix = ['sub-', options.subj.subjId, '_desc-'];
    for tt=1:length(tissuelabels)
        savestl(nmesh, emesh(emesh(:,5)==tt,1:4), fullfile(headmodelDir, [filePrefix, 'headmodel', num2str(side), '_label-', tissuelabels{tt}, '.stl']), tissuelabels{tt});
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
