function varargout=ea_genvat_simbio_ws(varargin)
% This function generates a volume of activated tissue around for each
% electrode.
% Usage: VAT=ea_genvat(coords_mm,stimparams,options).
% ? stimparams is a struct variable with fields U (8*1 with voltage
% entries) and Im (8*1 with Impedance measurements).
% 
% This function only touches the .VAT entry of stimparams struct of the
% given side.

if nargin==4
acoords=varargin{1};
stimparams=varargin{2};
side=varargin{3};
options=varargin{4};
elseif nargin==1

if ischar(varargin{1}) % return name of method.
    varargout{1}='Simbio_ws';
    return
end
end

%load('empirical_testdata'); % will produce data taken from lead dbs: 'coords','stimparams','side','options'

for side=options.sides
coords=acoords{side};
options.earoot=[fileparts(which('lead')),filesep];

%% some preprocessing to establish the lead trajectory
traj=[coords(4,:)+(coords(4,:)-coords(2,:));coords(1,:)+(coords(1,:)-coords(3,:))];
for d=1:3
itraj(:,d)=linspace(traj(1,d),traj(2,d));
end

%% convert trajectory mm2vox
V=spm_vol([options.earoot,'atlases',filesep,options.atlasset,filesep,'gm_mask.nii']);
trajmm=[itraj,ones(length(itraj),1)];
trajvox=V.mat\trajmm';
trajvox=trajvox(1:3,:)';


%% we will now produce a cubic headmodel that is aligned around the electrode using lead dbs:
[cimat,~,mat]=ea_sample_cuboid(trajvox,options,[options.earoot,'atlases',filesep,options.atlasset,filesep,'gm_mask.nii'],0,50,101);
mat=mat';
Vexp=ea_synth_nii([options.root,options.patientname,filesep,'tmp.nii'],mat,[2,0],cimat);
spm_write_vol(Vexp,cimat);


%% get electrodes handles:
resultfig=getappdata(gcf,'resultfig');
el_render=getappdata(resultfig,'el_render');
thiselhandle=el_render(1).el_render{side};
% fields: 1: trajectory body; 2: trajectory bottom; 3: trajectory top
% next three: contact one, etc.
% next three: contact spacing one, etc.
% last: tip

% establish coordinate grid:
nii=ea_load_nii([options.root,options.patientname,filesep,'tmp.nii']);
[xx,yy,zz]=ind2sub(size(nii.img),1:numel(nii.img));
XYZvx=[xx;yy;zz;ones(1,length(xx))];
XYZmm=Vexp.mat*XYZvx;
clear XYZvx

cnt=1;
ea_dispercent(0,'Exporting electrode components');
Xcon=nii.img; Xcon(:)=0; % initialize image for all contacts
Xins=nii.img; Xins(:)=0; % initialize image for all insulated electrode parts

for comp=1:options.elspec.numel*2+1
    ea_dispercent(comp/(options.elspec.numel*2+1));
    try % shaft and contacts, here three surface components
        cyl=thiselhandle(cnt); top=thiselhandle(cnt+1); bottom=thiselhandle(cnt+2);
        cyl = surf2patch(cyl,'triangles');
        cyl.faces=[cyl.faces;top.Faces+length(cyl.vertices);bottom.Faces+length(cyl.vertices)+length(top.Vertices)];
        cyl.vertices=[cyl.vertices;top.Vertices;bottom.Vertices];
        cnt=cnt+3;
    catch % tip of the electrode, here only one surface component..
        cyl=thiselhandle(cnt);
        cyl = surf2patch(cyl,'triangles'); 
    end
in=ea_intriangulation(cyl.vertices,cyl.faces,XYZmm(1:3,:)');
Xt=nii.img;
Xt(:)=0; Xt(in)=1;
if comp>1 && comp<options.elspec.numel+2 % these are the CONTACTS
Xcon=Xcon+Xt;
else % these are the insulated shaft, tip and spacings..
    Xins=Xins+Xt;
end
end

ea_dispercent(1,'end');

% for now, dipole is always placed at contact no. 2..
dpvx=Vexp.mat\[coords(2,:),1]';
dpvx=dpvx(1:3)';

%% read in gm data and convert to segmented mri
smri=ft_read_mri(Vexp.fname,'dataformat','nifti_spm');
smri.unit='mm';
smri=rmfield(smri,'anatomy');
smri.gray=logical(cimat);
smri.white=~smri.gray;
smri.contacts=logical(Xcon);
smri.insulation=logical(Xins);
smri.gray(smri.contacts)=0; smri.gray(smri.insulation)=0;
smri.white(smri.contacts)=0; smri.white(smri.insulation)=0;

%% create the mesh using fieldtrip:
cfg        = [];
cfg.tissue      = {'gray','white','contacts','insulation'};
cfg.method = 'hexahedral';
mesh       = ft_prepare_mesh(cfg,smri);
mesh = ft_transform_geometry(inv(smri.transform), mesh);



% plot gray matter:
% pmesh=mesh;
% pmesh.hex=pmesh.hex(pmesh.tissue==1,:);
% figure, ft_plot_mesh(pmesh,'vertexcolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.2,'edgealpha',0.1);

%% volume conductor
vol=ft_headmodel_simbio(mesh,'conductivity',[0.33 0.14 0.999 0.001]);
%% 




clear XYZvx
% create a grid of sensors around the dipole..
cnt=1;
dist=2;
for xx=-dist:dist
    for yy=-dist:dist
        for zz=-dist:dist
         XYZvx(cnt,:)=dpvx+[xx,yy,zz];
        cnt=cnt+1;
        end
    end
end



sens.elecpos=XYZvx(:,1:3);
% sensors need labels, so i label them from 1 to numel.
sens.label=arrayfun(@num2str,1:size(sens.elecpos,1),'UniformOutput',0);
% the following needs to be added to prevent an error.
sens.unit='vox';
vol.unit='vox';

vol=ft_prepare_vol_sens(vol,sens);



%dpmm=[coords(2,:),1];
%dpvx=Vexp.mat\dpmm';
%dpvx=dpvx(1:3,:)';


lf=leadfield_simbio(dpvx,vol);

keyboard

% plot lead-field:
figure
plot3(dpvx(1),dpvx(2),dpvx(3),'r*');
hold on


quiver3(sens.elecpos(:,1),sens.elecpos(:,2),sens.elecpos(:,3),lf(:,1),lf(:,2),lf(:,3));
end



function Vexp=ea_synth_nii(fname,mat,dt,Y)
Vexp.mat=mat;
Vexp.fname=fname;
Vexp.dt=dt;
Vexp.dim=size(Y);
Vexp.n=[1 1];
Vexp.descrip='lead dbs - leadfield';
Vexp.private=nifti;
Vexp.private.mat=Vexp.mat;
Vexp.private.mat0=Vexp.mat;
Vexp.private.mat_intent='Aligned';
Vexp.private.mat0_intent='Aligned';
Vexp.private.descrip='lead dbs - leadfield';
%Vexp.private.dat=Y;





