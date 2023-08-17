function ea_segment_electrode(~,~,options,resultfig,onoff)

directory=[options.root,options.patientname,filesep];
if options.native
    switch options.subj.postopModality
        case 'MRI'
            elnii = options.subj.coreg.anat.postop.ax_MRI;
        case 'CT'
            elnii = options.subj.coreg.anat.postop.CT;
    end
    elssubf='native';
else
    switch options.subj.postopModality
        case 'MRI'
            elnii=options.subj.norm.anat.postop.ax_MRI;
        case 'CT'
            elnii=options.subj.norm.anat.postop.CT;
    end
    elssubf='template';
end

switch options.subj.postopModality
    case 'MRI'
        tval=-50;
    case 'CT'
        tval=2500;
end

switch onoff
    case 'on'
        % check if has been visualized before
        elseg=getappdata(resultfig,'elseg');
        if ~isempty(elseg) % && 0
            elseg.Visible='on';
        else
            % check if segmentation exists
            nii=ea_load_nii(elnii);

            if options.native
                if exist(options.subj.brainshift.transform.scrf,'file') % apply brainshift correction to files on the fly.
                    scrf=load(options.subj.brainshift.transform.scrf);
                    nii.mat=scrf.mat*nii.mat;
                end
            end

            if strcmp(options.subj.postopModality, 'MRI') % not yet implemented for MR.
                nii.img=-nii.img;
            end

            fv=isosurface(permute(nii.img,[2,1,3]),tval);
            fvc=isocaps(permute(nii.img,[2,1,3]),tval);
            fv.faces=[fv.faces;fvc.faces+size(fv.vertices,1)];
            fv.vertices=[fv.vertices;fvc.vertices];

            % fv=isosurface(nii.img,2500); % only CT support for now
            % figure, patch('faces',fv.faces,'vertices',fv.vertices,'Edgecolor','none')
            fv.vertices=[fv.vertices,ones(size(fv.vertices,1),1)];
            fv.vertices=fv.vertices*nii.mat';

            %check and apply brainshift correction:
            if exist([directory,'scrf',filesep,'scrf_converted.mat'],'file') && options.native
                d=load([directory,'scrf',filesep,'scrf_converted.mat']);
                bsmat=d.mat;
                fv.vertices=fv.vertices*bsmat';
            end

            fv.vertices=fv.vertices(:,1:3);
            % % save fv
            % if ~exist([directory,'electrode_auto_segment',filesep,elssubf,filesep],'dir')
            %     mkdir([directory,'electrode_auto_segment',filesep,elssubf,filesep]);
            % end

            fv=refinepatch(fv);
            fv=ea_smoothpatch(fv,0,10);

            elseg=patch('faces',fv.faces,'vertices',fv.vertices,'Edgecolor','none','FaceColor','interp','FaceVertexCData',repmat([0.9,0.9,0.8],size(fv.vertices,1),1));

            setappdata(resultfig,'elseg',elseg);
        end
    case 'off'
        elseg=getappdata(resultfig,'elseg');
        if ~isempty(elseg)
            elseg.Visible='off';
        end
end

function [FV2]=refinepatch(FV)
% This function "refinepatch" refines a triangular mesh with
% a spline interpolated 4-split method.
%
%   [FV2] = refinepatch(FV,options)
%
% inputs,
%   FV : Structure containing a Patch, with
%        FV.vertices the mesh vertices
%        FV.face the mesh faces (triangles), rows with each 3 vertex indices
% outputs,
%   FV2 : Structure Containing the refined patch
%
%
% Reference:
%  The spline interpolation of the face edges is done by the
%  Opposite Edge Method, described in: "Construction of Smooth Curves
%  and Surfaces from Polyhedral Models" by Leon A. Shirman
%
% How it works:
%  The tangents (normals) and velocity on the edge points of all edges
%  are calculated. Which are  later used for b-spline interpolation when
%  splitting the edges.
%
%  A tangent on an 3D line or edge is under defined and can rotate along
%  the line, thus an (virtual) opposite vertex is used to fix the tangent and
%  make it more like a surface normal.
%
%  B-spline interpolate a half way vertices between all existing vertices
%  using the velocity and tangent from the edgepoints. After splitting a
%  new facelist is constructed
%
% Speed:
%  Compile the c-functions for more speed with:
%   mex vertex_neighbours_double.c -v;
%   mex edge_tangents_double.c -v;
%
% Example:
%
% X=[-0.5000;  0.5000;  0.0000;  0.0000];
% Y=[-0.2887; -0.2887;  0.5774;  0.0000];
% Z=[ 0.0000;  0.0000;  0.0000;  0.8165];
% FV.vertices=[X Y Z];
%
% FV.faces=[2 3 4; 4 3 1; 1 2 4; 3 2 1];
%
% figure, set(gcf, 'Renderer', 'opengl'); axis equal;
% for i=1:4
%   patch(FV,'facecolor',[1 0 0]);
%   pause(2);
%   [FV]=refinepatch(FV);
% end
%
% Function is written by D.Kroon University of Twente (February 2010)

% Get the neighbour vertices of each vertice from the face list.
Ne=vertex_neighbours(FV);

% Calculate the tangents (normals) and velocity of all edges. Which is
% later used for b-spline interpolation and split of the edges
%
% A tangent on an 3D line or edge is under defined and can rotate along
% the line, thus an (virtual) opposite vertex is used to fix the tangent and
% make it more like a surface normal.
V=FV.vertices; F=FV.faces;
[ET_table,EV_table,ETV_index]=edge_tangents(V,Ne);

% B-spline interpolate a half way vertices between all existing vertices
% using the velocity and tangent from above
[V,HT_index, HT_values]=make_halfway_vertices(EV_table,ET_table,ETV_index,V,Ne);

% Make new facelist
Fnew=makenewfacelist(F,HT_index,HT_values);

FV2.vertices=V;
FV2.faces=Fnew;

function Ne=vertex_neighbours(FV)
% This function VERTEX_NEIGHBOURS will search in a face list for all
% the neigbours of each vertex.
%
% Ne=vertex_neighbours(FV)
%

Ne=vertex_neighbours_double(FV.faces(:,1),FV.faces(:,2),FV.faces(:,3),FV.vertices(:,1),FV.vertices(:,2),FV.vertices(:,3));


function Ne=vertex_neighbours_double(Fa,Fb,Fc,Vx,Vy,Vz)

F=[Fa Fb Fc];
V=[Vx Vy Vz];

% Neighbourh cell array
Ne=cell(1,size(V,1));

% Loop through all faces
for i=1:length(F)
    % Add the neighbors of each vertice of a face
    % to his neighbors list.
    Ne{F(i,1)}=[Ne{F(i,1)} [F(i,2) F(i,3)]];
    Ne{F(i,2)}=[Ne{F(i,2)} [F(i,3) F(i,1)]];
    Ne{F(i,3)}=[Ne{F(i,3)} [F(i,1) F(i,2)]];
end

% Loop through all neighbor arrays and sort them (Rotation same as faces)
for i=1:size(V,1)

    Pneighf=Ne{i};
    if(isempty(Pneighf))
        Pneig=[];
    else
        start=1;
        for index1=1:2:length(Pneighf)
            found=false;
            for index2=2:2:length(Pneighf),
                if(Pneighf(index1)==Pneighf(index2))
                    found=true; break
                end
            end
            if(~found)
                start=index1; break
            end
        end
        Pneig=[];
        Pneig(1)=Pneighf(start);
        Pneig(2)=Pneighf(start+1);

        % Add the neighbours with respect to original rotation
        for j=2+double(found):(length(Pneighf)/2)
            found = false;
            for index=1:2:length(Pneighf),
                if(Pneighf(index)==Pneig(end))
                    if(sum(Pneig==Pneighf(index+1))==0)
                        found =true;
                        Pneig=[Pneig Pneighf(index+1)];
                    end
                end
            end
            if(~found) % This only happens with weird edge vertices
                for index=1:2:length(Pneighf),
                    if(sum(Pneig==Pneighf(index))==0)
                        Pneig=[Pneig Pneighf(index)];
                        if(sum(Pneig==Pneighf(index+1))==0)
                            Pneig=[Pneig Pneighf(index+1)];
                        end
                    end
                end
            end
        end
        % Add forgotten neigbours
        if(length(Pneig)<length(Pneighf))
            for j=1:length(Pneighf)
                if(sum(Pneig==Pneighf(j))==0)
                    Pneig=[Pneig Pneighf(j)];
                end
            end
        end
    end
    Ne{i}=Pneig;
end


function [ET_table,EV_table,ETV_index]=edge_tangents(V,Ne)
[ET_table,EV_table,ETV_index]=edge_tangents_double(double(V),Ne);

function [ET_table,EV_table,ETV_index]=edge_tangents_double(V,Ne)
% Edge tangents table
ET_table=zeros(size(V,1)*4,3);
% Edge velocity table
EV_table=zeros(size(V,1)*4,1);
% Edge tangents/velocity index for tables
ETV_index=zeros(size(V,1)*4,2);
ETV_num=0;

% Calculate the tangents and velocity for each edge
Pn=zeros(length(Ne),3); Pnop=zeros(length(Ne),3);
for i=1:size(V,1)
    P=V(i,:);
    Pneig=Ne{i};

    % Find the opposite vertex of each neigbourh vertex.
    % incase of odd number of neigbourhs interpolate the opposite neigbourh
    if(mod(length(Pneig),2)==0)
        for k=1:length(Pneig)
            neg=k+length(Pneig)/2; if(neg>length(Pneig)), neg=neg-length(Pneig); end
            Pn(k,:) = V(Pneig(k),:); Pnop(k,:) = V(Pneig(neg),:);
        end
    else
        for k=1:length(Pneig)
            neg=k+length(Pneig)/2; neg1=floor(neg); neg2=ceil(neg);
            if(neg1>length(Pneig)), neg1=neg1-length(Pneig); end
            if(neg2>length(Pneig)), neg2=neg2-length(Pneig); end
            Pn(k,:) = V(Pneig(k),:); Pnop(k,:) = (V(Pneig(neg1),:)+V(Pneig(neg2),:))/2;
        end
    end

    for j=1:length(Pneig);
        % Calculate length edges of face
        Ec=sqrt(sum((Pn(j,:)-P).^2))+1e-14;
        Eb=sqrt(sum((Pnop(j,:)-P).^2))+1e-14;
        Ea=sqrt(sum((Pn(j,:)-Pnop(j,:)).^2))+1e-14;

        % Calculate face surface area
        s = ((Ea+Eb+Ec)/2);
        h = (2/Ea)*sqrt(s*(s-Ea)*(s-Eb)*(s-Ec))+1e-14;
        x = (Ea^2-Eb^2+Ec^2)/(2*Ea);

        % 2D triangle coordinates
        % corx(1)=0;    cory(1)=0;
        % corx(2)=x;    cory(2)=h;
        % corx(3)=Ea;   cory(3)=0;
        % corx(4)=0;    cory(4)=0;

        % Calculate tangent of 2D triangle
        Np=[-h x]; Np=Np/(sqrt(sum(Np.^2))+1e-14);
        Ns=[h Ea-x]; Ns=Ns/(sqrt(sum(Ns.^2))+1e-14);
        Nb=Np+Ns;
        Tb=[Nb(2) -Nb(1)];

        % Back to 3D coordinates
        Pm=(Pn(j,:)*x+Pnop(j,:)*(Ea-x))/Ea;
        X3=(Pn(j,:)-Pnop(j,:))/Ea;
        Y3=(P-Pm)/h;

        % 2D tangent to 3D tangent
        Tb3D=(X3*Tb(1)+Y3*Tb(2));  Tb3D=Tb3D/(sqrt(sum(Tb3D.^2))+1e-14);

        % Edge Velocity
        Vv=0.5*(Ec+0.5*Ea);

        ETV_num=ETV_num+1;
        ETV_index(ETV_num,:)=[i Pneig(j)];
        ET_table(ETV_num,:)= Tb3D;
        EV_table(ETV_num)=Vv;
    end
end


function [Vout,HT_index, HT_values]=make_halfway_vertices(EV_table,ET_table,ETV_index,V,Ne)
% ETV_index_sparse=sparse(size(V,1),size(V,1));
%
% for i=1:size(ETV_index,1)
% 	if(ETV_index(i,1)>0)
% 		ETV_index_sparse(ETV_index(i,1),ETV_index(i,2))=i;
% 	else
% 		ET_table=ET_table(1:(i-1),:);
% 		EV_table=EV_table(1:(i-1));
% 		break
% 	end
% end


% Table to cell array
ETV_index_vall=cell(length(V),1);
for i=1:size(ETV_index,1)
    if(ETV_index(i,1)>0)
        ETV_index_vall{ETV_index(i,1)}=[ETV_index_vall{ETV_index(i,1)} i];
    else
        ET_table=ET_table(1:(i-1),:);
        EV_table=EV_table(1:(i-1));
        break;
    end
end

HT_index=cell(length(V),1);
HT_values=cell(length(V),1);

% Make output V
Vout=zeros(size(V,1)*4,3);
Vout(1:size(V,1),:)=V;
Vindex=size(V,1);

for i=1:length(V)
    Pneig=Ne{i};
    for j=1:length(Pneig);
        % Get the tangent and velocity of the edge P -> Pneig
        index=Ne{i}; vals=ETV_index_vall{i}; select1=vals(index==Pneig(j));
        Va=EV_table( select1); Ea=ET_table( select1,:);
        % Get the tangent and velocity of the edge Pneig -> P
        index=Ne{Pneig(j)}; vals=ETV_index_vall{Pneig(j)}; select2=vals(index==i);
        Vb=EV_table(select2); Eb=ET_table(select2,:);

        % The four points describing the spline
        P0=V(i,:);
        P3=V(Pneig(j),:);
        P1=P0+Ea*Va/3;
        P2=P3+Eb*Vb/3;

        % Spline used to calculated the xyz coordinate of the middle of each edge;
        c = 3*(P1 - P0);
        b = 3*(P2 - P1) - c;
        a = P3 - P0 - c - b;

        halfwayp = a*0.125 + b*0.250 + c*0.500 + P0;

        % Save the edge middle point
        if(sum(HT_index{i}==Pneig(j))==0)
            Vindex=Vindex+1;
            Vout(Vindex,:)=halfwayp;
            HT_index {i}= [HT_index{i} Pneig(j)];
            HT_values{i}=[HT_values{i} Vindex];
            HT_index {Pneig(j)}=[HT_index{ Pneig(j)} i];
            HT_values{Pneig(j)}=[HT_values{Pneig(j)} Vindex];
        end
    end
end
Vout=Vout(1:Vindex,:);

function Fnew=makenewfacelist(F,HT_index, HT_values)
% Combine the edge middle points and old vertex points to faces.
% (4 Split method)
Fnew=zeros(length(F)*4,3);
for i=0:length(F)-1,
    vert1=F(i+1,1);
    vert2=F(i+1,2);
    vert3=F(i+1,3);

    index=HT_index{vert1}; vals=HT_values{vert1};
    verta= vals(index==vert2);
    index=HT_index{vert2}; vals=HT_values{vert2};
    vertb= vals(index==vert3);
    index=HT_index{vert3}; vals=HT_values{vert3};
    vertc= vals(index==vert1);

    Fnew(i*4+1,:)=[vert1 verta vertc];
    Fnew(i*4+2,:)=[verta vert2 vertb];
    Fnew(i*4+3,:)=[vertc vertb vert3];
    Fnew(i*4+4,:)=[verta vertb vertc];
end
