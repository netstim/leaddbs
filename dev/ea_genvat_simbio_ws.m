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


usediffusion=0; % set to 1 to incorporate diffusion signal (for now only possible using the mesoFT tracker).

%load('empirical_testdata'); % will produce data taken from lead dbs: 'coords','stimparams','side','options'

for side=options.sides
coords=acoords{side};
options.earoot=[ea_getearoot];

%% some preprocessing to establish the lead trajectory
traj=[coords(4,:)+(coords(4,:)-coords(2,:));coords(1,:)+(coords(1,:)-coords(3,:))];
for d=1:3
itraj(:,d)=linspace(traj(1,d),traj(2,d));
end

%% convert trajectory mm2vox
V=spm_vol([ea_space(options,'atlases'),options.atlasset,filesep,'gm_mask.nii']);
trajmm=[itraj,ones(length(itraj),1)];
trajvox=V.mat\trajmm';
trajvox=trajvox(1:3,:)';


%% we will now produce a cubic headmodel that is aligned around the electrode using lead dbs:
[cimat,~,mat]=ea_sample_cuboid(trajvox,options,[ea_space(options,'atlases'),options.atlasset,filesep,'gm_mask.nii'],0,25,51); % this will result in ~10x10x10 mm.
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
    % this following method takes quite some time... even more importantly,
    % the info will be transfered from mesh to volume and lateron back to
    % mesh again. For now, this is still the most convenient method.
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

% set up dipole


% if length(find(stimparams.U))>1 % bipolar
% elseif length(find(stimparams.U)==1 % monopolar
% else % no active contact!
%     ea_error('No active stimulation contact selected.');
% end

% define dipole coordinates:
    dpvx=coords(find(stimparams.U),:);

dpvx=Vexp.mat\[dpvx,ones(size(dpvx,1),1)]';
dpvx=dpvx(1:3,:)';

%% read in gm data and convert to segmented mri
% construct ft-like anatomy structure based on SPM nifti info.
smri.dim=Vexp.dim;
smri.hdr=Vexp;
smri.transform=Vexp.mat;
smri.unit='mm';
smri.gray=logical(cimat); % gm portion of the file
smri.white=~smri.gray; % first, set everything that is not gm to wm
smri.contacts=logical(Xcon); % include electrode contacts to the model
smri.insulation=logical(Xins); % include insulation to the model
smri.gray(smri.contacts)=0; smri.gray(smri.insulation)=0; % remove contact and insulation portions from the gm..
smri.white(smri.contacts)=0; smri.white(smri.insulation)=0; % .. and white matter portions.

% export a nifti version of the headmodel just for control reasons and as
% an export for other use-cases:
X=single(smri.gray);
X=X+2*smri.white;
X=X+3*smri.contacts;
X=X+4*smri.insulation;
mkdir([options.root,options.patientname,filesep,'headmodel']);
Vex=Vexp; Vex.fname=[options.root,options.patientname,filesep,'headmodel',filesep,'structural.nii'];
spm_write_vol(Vex,X);
delete([options.root,options.patientname,filesep,'tmp.nii']);
clear Vex X

%% generate diffusion signal:
if usediffusion
    disp('Loading FTR...');
    load([options.earoot,'dev',filesep,'bTensor']) % b-Tensor of a simple 6fold diffusion series
    ftr=load([options.root,options.patientname,filesep,options.prefs.FTR_normalized]);
    disp('Done. Estimating diffusion signal based on fibertracts...');
    signal=ea_ftr2Sigmaps(ftr,ten);
    disp('Done. Calculating Tensors...');
    reftemplate=[ea_space(options,'dartel'),'dartelmni_1.nii,2'];
    Vsig=spm_vol(reftemplate);
    for i=1:size(signal,4);
        Vsig.fname=[options.root,options.patientname,filesep,'headmodel',filesep,'dti_',num2str(i),'.nii'];
        spm_write_vol(Vsig,squeeze(signal(:,:,:,i)));
    end

end


%% create the mesh using fieldtrip:
cfg        = [];
cfg.tissue      = {'gray','white','contacts','insulation'};
cfg.method = 'hexahedral';
mesh       = ea_ft_prepare_mesh(cfg,smri); % need to incorporate this function and all dependencies into lead-dbs
mesh = ea_ft_transform_geometry(inv(smri.transform), mesh);



% plot gray matter:
% pmesh=mesh;
% pmesh.hex=pmesh.hex(pmesh.tissue==1,:);
% figure, ft_plot_mesh(pmesh,'vertexcolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.2,'edgealpha',0.1);

%% volume conductor
vol=ea_ft_headmodel_simbio(mesh,'conductivity',[0.33 0.14 0.999 0.001]); % need to incorporate this function and all dependencies into lead-dbs
%%




clear XYZvx
% create a grid of sensors around the dipole..
cnt=1;
dist=7; % how far to leave the dipole in each direction (in voxels)
swidth=2; % step width (in voxels)
for xx=-dist:swidth:dist
    for yy=-dist:swidth:dist
        for zz=-dist:swidth:dist
         XYZvx(cnt,:)=mean(dpvx,1)+[xx,yy,zz];
        cnt=cnt+1;
        end
    end
end
sens.elecpos=XYZvx(:,1:3);
% sensors need labels, so we label them from 1 to numel.
sens.label=arrayfun(@num2str,1:size(sens.elecpos,1),'UniformOutput',0);
sens.unit='vox';
vol.unit='vox';

vol=ea_ft_prepare_vol_sens(vol,sens);
lf=ea_leadfield_simbio(dpvx,vol);



% plot lead-field:
figure
plot3(dpvx(:,1),dpvx(:,2),dpvx(:,3),'r*');
hold on


quiver3(sens.elecpos(:,1),sens.elecpos(:,2),sens.elecpos(:,3),lf(:,1),lf(:,2),lf(:,3));
end


%% begin fieldtrip/simbio functions:

function [warped] = ea_ft_warp_apply(M, input, method, tol)

% FT_WARP_APPLY performs a 3D linear or nonlinear transformation on the input
% coordinates, similar to those in AIR 3.08. You can find technical
% documentation on warping in general at http://bishopw.loni.ucla.edu/AIR3
%
% Use as
%   [warped] = ft_warp_apply(M, input, method, tol)
% where
%   M        vector or matrix with warping parameters
%   input    Nx3 matrix with coordinates
%   warped   Nx3 matrix with coordinates
%   method   string describing the warping method
%   tol      (optional) value determining the numerical precision of the
%             output, to deal with numerical round off imprecisions due to
%             the warping
%
% The methods 'nonlin0', 'nonlin2' ... 'nonlin5' specify a
% polynomial transformation. The size of the transformation matrix
% depends on the order of the warp
%   zeroth order :  1 parameter  per coordinate (translation)
%   first  order :  4 parameters per coordinate (total 12, affine)
%   second order : 10 parameters per coordinate
%   third  order : 20 parameters per coordinate
%   fourth order : 35 parameters per coordinate
%   fifth  order : 56 parameters per coordinate (total 168)
% The size of M should be 3xP, where P is the number of parameters
% per coordinate. Alternatively, you can specify the method to be
% 'nonlinear', where the order will be determined from the size of
% the matrix M.
%
% If the method 'homogeneous' is selected, the input matrix M should be
% a 4x4 homogenous transformation matrix.
%
% If the method 'sn2individual' or 'individual2sn' is selected, the input
% M should be a structure based on nonlinear (warping) normalisation parameters
% created by SPM8 for alignment between an individual structural MRI and the
% template MNI brain.  These options call private functions of the same name.
% M will have subfields like this:
%     Affine: [4x4 double]
%         Tr: [4-D double]
%         VF: [1x1 struct]
%         VG: [1x1 struct]
%      flags: [1x1 struct]
%
% If any other method is selected, it is assumed that it specifies
% the name of an auxiliary function that will, when given the input
% parameter vector M, return an 4x4 homogenous transformation
% matrix. Supplied functions in the warping toolbox are translate,
% rotate, scale, rigidbody, globalrescale, traditional, affine,
% perspective.
%
% See also FT_WARP_OPTIM, FT_WARP_ERROR

% Copyright (C) 2000-2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_warp_apply.m 10132 2015-01-27 16:08:29Z johzum $

if nargin<4
  tol = [];
end

if nargin<3 && all(size(M)==4)
  % no specific transformation mode has been selected
  % it looks like a homogenous transformation matrix
  method = 'homogeneous';
elseif nargin<3
  % the default method is 'nonlinear'
  method = 'nonlinear';
end

if size(input,2)==2
  % convert the input points from 2D to 3D representation
  input(:,3) = 0;
  input3d = false;
else
  input3d = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonlinear warping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp(method, {'nonlinear', 'nonlin0', 'nonlin1', 'nonlin2', 'nonlin3', 'nonlin4', 'nonlin5'}))
  x = input(:,1);
  y = input(:,2);
  z = input(:,3);
  s = size(M);

  if s(1)~=3
    error('invalid size of nonlinear transformation matrix');
  elseif strcmp(method, 'nonlin0') && s(2)~=1
    error('invalid size of nonlinear transformation matrix');
  elseif strcmp(method, 'nonlin1') && s(2)~=4
    error('invalid size of nonlinear transformation matrix');
  elseif strcmp(method, 'nonlin2') && s(2)~=10
    error('invalid size of nonlinear transformation matrix');
  elseif strcmp(method, 'nonlin3') && s(2)~=20
    error('invalid size of nonlinear transformation matrix');
  elseif strcmp(method, 'nonlin4') && s(2)~=35
    error('invalid size of nonlinear transformation matrix');
  elseif strcmp(method, 'nonlin5') && s(2)~=56
    error('invalid size of nonlinear transformation matrix');
  end

  if s(2)==1
    % this is a translation, which in a strict sense is not the 0th order nonlinear transformation
    xx = M(1,1) + x;
    yy = M(2,1) + y;
    zz = M(3,1) + z;
  elseif s(2)==4
    xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z;
    yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z;
    zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z;
  elseif s(2)==10
    xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z;
    yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z;
    zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z;
  elseif s(2)==20
    xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z + M(1,11)*x.*x.*x + M(1,12)*x.*x.*y + M(1,13)*x.*x.*z + M(1,14)*x.*y.*y + M(1,15)*x.*y.*z + M(1,16)*x.*z.*z + M(1,17)*y.*y.*y + M(1,18)*y.*y.*z + M(1,19)*y.*z.*z + M(1,20)*z.*z.*z;
    yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z + M(2,11)*x.*x.*x + M(2,12)*x.*x.*y + M(2,13)*x.*x.*z + M(2,14)*x.*y.*y + M(2,15)*x.*y.*z + M(2,16)*x.*z.*z + M(2,17)*y.*y.*y + M(2,18)*y.*y.*z + M(2,19)*y.*z.*z + M(2,20)*z.*z.*z;
    zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z + M(3,11)*x.*x.*x + M(3,12)*x.*x.*y + M(3,13)*x.*x.*z + M(3,14)*x.*y.*y + M(3,15)*x.*y.*z + M(3,16)*x.*z.*z + M(3,17)*y.*y.*y + M(3,18)*y.*y.*z + M(3,19)*y.*z.*z + M(3,20)*z.*z.*z;
  elseif s(2)==35
    xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z + M(1,11)*x.*x.*x + M(1,12)*x.*x.*y + M(1,13)*x.*x.*z + M(1,14)*x.*y.*y + M(1,15)*x.*y.*z + M(1,16)*x.*z.*z + M(1,17)*y.*y.*y + M(1,18)*y.*y.*z + M(1,19)*y.*z.*z + M(1,20)*z.*z.*z + M(1,21)*x.*x.*x.*x + M(1,22)*x.*x.*x.*y + M(1,23)*x.*x.*x.*z + M(1,24)*x.*x.*y.*y + M(1,25)*x.*x.*y.*z + M(1,26)*x.*x.*z.*z + M(1,27)*x.*y.*y.*y + M(1,28)*x.*y.*y.*z + M(1,29)*x.*y.*z.*z + M(1,30)*x.*z.*z.*z + M(1,31)*y.*y.*y.*y + M(1,32)*y.*y.*y.*z + M(1,33)*y.*y.*z.*z + M(1,34)*y.*z.*z.*z + M(1,35)*z.*z.*z.*z;
    yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z + M(2,11)*x.*x.*x + M(2,12)*x.*x.*y + M(2,13)*x.*x.*z + M(2,14)*x.*y.*y + M(2,15)*x.*y.*z + M(2,16)*x.*z.*z + M(2,17)*y.*y.*y + M(2,18)*y.*y.*z + M(2,19)*y.*z.*z + M(2,20)*z.*z.*z + M(2,21)*x.*x.*x.*x + M(2,22)*x.*x.*x.*y + M(2,23)*x.*x.*x.*z + M(2,24)*x.*x.*y.*y + M(2,25)*x.*x.*y.*z + M(2,26)*x.*x.*z.*z + M(2,27)*x.*y.*y.*y + M(2,28)*x.*y.*y.*z + M(2,29)*x.*y.*z.*z + M(2,30)*x.*z.*z.*z + M(2,31)*y.*y.*y.*y + M(2,32)*y.*y.*y.*z + M(2,33)*y.*y.*z.*z + M(2,34)*y.*z.*z.*z + M(2,35)*z.*z.*z.*z;
    zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z + M(3,11)*x.*x.*x + M(3,12)*x.*x.*y + M(3,13)*x.*x.*z + M(3,14)*x.*y.*y + M(3,15)*x.*y.*z + M(3,16)*x.*z.*z + M(3,17)*y.*y.*y + M(3,18)*y.*y.*z + M(3,19)*y.*z.*z + M(3,20)*z.*z.*z + M(3,21)*x.*x.*x.*x + M(3,22)*x.*x.*x.*y + M(3,23)*x.*x.*x.*z + M(3,24)*x.*x.*y.*y + M(3,25)*x.*x.*y.*z + M(3,26)*x.*x.*z.*z + M(3,27)*x.*y.*y.*y + M(3,28)*x.*y.*y.*z + M(3,29)*x.*y.*z.*z + M(3,30)*x.*z.*z.*z + M(3,31)*y.*y.*y.*y + M(3,32)*y.*y.*y.*z + M(3,33)*y.*y.*z.*z + M(3,34)*y.*z.*z.*z + M(3,35)*z.*z.*z.*z;
  elseif s(2)==56
    xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z + M(1,11)*x.*x.*x + M(1,12)*x.*x.*y + M(1,13)*x.*x.*z + M(1,14)*x.*y.*y + M(1,15)*x.*y.*z + M(1,16)*x.*z.*z + M(1,17)*y.*y.*y + M(1,18)*y.*y.*z + M(1,19)*y.*z.*z + M(1,20)*z.*z.*z + M(1,21)*x.*x.*x.*x + M(1,22)*x.*x.*x.*y + M(1,23)*x.*x.*x.*z + M(1,24)*x.*x.*y.*y + M(1,25)*x.*x.*y.*z + M(1,26)*x.*x.*z.*z + M(1,27)*x.*y.*y.*y + M(1,28)*x.*y.*y.*z + M(1,29)*x.*y.*z.*z + M(1,30)*x.*z.*z.*z + M(1,31)*y.*y.*y.*y + M(1,32)*y.*y.*y.*z + M(1,33)*y.*y.*z.*z + M(1,34)*y.*z.*z.*z + M(1,35)*z.*z.*z.*z + M(1,36)*x.*x.*x.*x.*x + M(1,37)*x.*x.*x.*x.*y + M(1,38)*x.*x.*x.*x.*z + M(1,39)*x.*x.*x.*y.*y + M(1,40)*x.*x.*x.*y.*z + M(1,41)*x.*x.*x.*z.*z + M(1,42)*x.*x.*y.*y.*y + M(1,43)*x.*x.*y.*y.*z + M(1,44)*x.*x.*y.*z.*z + M(1,45)*x.*x.*z.*z.*z + M(1,46)*x.*y.*y.*y.*y + M(1,47)*x.*y.*y.*y.*z + M(1,48)*x.*y.*y.*z.*z + M(1,49)*x.*y.*z.*z.*z + M(1,50)*x.*z.*z.*z.*z + M(1,51)*y.*y.*y.*y.*y + M(1,52)*y.*y.*y.*y.*z + M(1,53)*y.*y.*y.*z.*z + M(1,54)*y.*y.*z.*z.*z + M(1,55)*y.*z.*z.*z.*z + M(1,56)*z.*z.*z.*z.*z;
    yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z + M(2,11)*x.*x.*x + M(2,12)*x.*x.*y + M(2,13)*x.*x.*z + M(2,14)*x.*y.*y + M(2,15)*x.*y.*z + M(2,16)*x.*z.*z + M(2,17)*y.*y.*y + M(2,18)*y.*y.*z + M(2,19)*y.*z.*z + M(2,20)*z.*z.*z + M(2,21)*x.*x.*x.*x + M(2,22)*x.*x.*x.*y + M(2,23)*x.*x.*x.*z + M(2,24)*x.*x.*y.*y + M(2,25)*x.*x.*y.*z + M(2,26)*x.*x.*z.*z + M(2,27)*x.*y.*y.*y + M(2,28)*x.*y.*y.*z + M(2,29)*x.*y.*z.*z + M(2,30)*x.*z.*z.*z + M(2,31)*y.*y.*y.*y + M(2,32)*y.*y.*y.*z + M(2,33)*y.*y.*z.*z + M(2,34)*y.*z.*z.*z + M(2,35)*z.*z.*z.*z + M(2,36)*x.*x.*x.*x.*x + M(2,37)*x.*x.*x.*x.*y + M(2,38)*x.*x.*x.*x.*z + M(2,39)*x.*x.*x.*y.*y + M(2,40)*x.*x.*x.*y.*z + M(2,41)*x.*x.*x.*z.*z + M(2,42)*x.*x.*y.*y.*y + M(2,43)*x.*x.*y.*y.*z + M(2,44)*x.*x.*y.*z.*z + M(2,45)*x.*x.*z.*z.*z + M(2,46)*x.*y.*y.*y.*y + M(2,47)*x.*y.*y.*y.*z + M(2,48)*x.*y.*y.*z.*z + M(2,49)*x.*y.*z.*z.*z + M(2,50)*x.*z.*z.*z.*z + M(2,51)*y.*y.*y.*y.*y + M(2,52)*y.*y.*y.*y.*z + M(2,53)*y.*y.*y.*z.*z + M(2,54)*y.*y.*z.*z.*z + M(2,55)*y.*z.*z.*z.*z + M(2,56)*z.*z.*z.*z.*z;
    zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z + M(3,11)*x.*x.*x + M(3,12)*x.*x.*y + M(3,13)*x.*x.*z + M(3,14)*x.*y.*y + M(3,15)*x.*y.*z + M(3,16)*x.*z.*z + M(3,17)*y.*y.*y + M(3,18)*y.*y.*z + M(3,19)*y.*z.*z + M(3,20)*z.*z.*z + M(3,21)*x.*x.*x.*x + M(3,22)*x.*x.*x.*y + M(3,23)*x.*x.*x.*z + M(3,24)*x.*x.*y.*y + M(3,25)*x.*x.*y.*z + M(3,26)*x.*x.*z.*z + M(3,27)*x.*y.*y.*y + M(3,28)*x.*y.*y.*z + M(3,29)*x.*y.*z.*z + M(3,30)*x.*z.*z.*z + M(3,31)*y.*y.*y.*y + M(3,32)*y.*y.*y.*z + M(3,33)*y.*y.*z.*z + M(3,34)*y.*z.*z.*z + M(3,35)*z.*z.*z.*z + M(3,36)*x.*x.*x.*x.*x + M(3,37)*x.*x.*x.*x.*y + M(3,38)*x.*x.*x.*x.*z + M(3,39)*x.*x.*x.*y.*y + M(3,40)*x.*x.*x.*y.*z + M(3,41)*x.*x.*x.*z.*z + M(3,42)*x.*x.*y.*y.*y + M(3,43)*x.*x.*y.*y.*z + M(3,44)*x.*x.*y.*z.*z + M(3,45)*x.*x.*z.*z.*z + M(3,46)*x.*y.*y.*y.*y + M(3,47)*x.*y.*y.*y.*z + M(3,48)*x.*y.*y.*z.*z + M(3,49)*x.*y.*z.*z.*z + M(3,50)*x.*z.*z.*z.*z + M(3,51)*y.*y.*y.*y.*y + M(3,52)*y.*y.*y.*y.*z + M(3,53)*y.*y.*y.*z.*z + M(3,54)*y.*y.*z.*z.*z + M(3,55)*y.*z.*z.*z.*z + M(3,56)*z.*z.*z.*z.*z;
  else
    error('invalid size of nonlinear transformation matrix');
  end

  warped = [xx yy zz];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % linear warping using homogenous coordinate transformation matrix
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(method, 'homogenous') || strcmp(method, 'homogeneous')
  if all(size(M)==3)
    % convert the 3x3 homogenous transformation matrix (corresponding with 2D)
    % into a 4x4 homogenous transformation matrix (corresponding with 3D)
    M = [
      M(1,1) M(1,2)  0  M(1,3)
      M(2,1) M(2,2)  0  M(2,3)
      0      0       0  0
      M(3,1) M(3,2)  0  M(3,3)
      ];
  end

  %warped = M * [input'; ones(1, size(input, 1))];
  %warped = warped(1:3,:)';

  % below achieves the same as lines 154-155
  warped = [input ones(size(input, 1),1)]*M(1:3,:)';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % using external function that returns a homogeneous transformation matrix
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif exist(method, 'file') && ~isa(M, 'struct')
  % get the homogenous transformation matrix
  H = feval(method, M);
  warped = ea_ft_warp_apply(H, input, 'homogeneous');

elseif strcmp(method, 'sn2individual') && isa(M, 'struct')
  % use SPM structure with parameters for an inverse warp
  % from normalized space to individual, can be non-linear
  warped = sn2individual(M, input);

elseif strcmp(method, 'individual2sn') && isa(M, 'struct')
  % use SPM structure with parameters for a warp from
  % individual space to normalized space, can be non-linear
  %error('individual2sn is not yet implemented');
  warped = individual2sn(M, input);
else
  error('unrecognized transformation method');
end

if ~input3d
  % convert from 3D back to 2D representation
  warped = warped(:,1:2);
end

if ~isempty(tol)
  if tol>0
    warped = fix(warped./tol)*tol;
  end
end

function [vol, sens] = ea_ft_prepare_vol_sens(vol, sens, varargin)

% FT_PREPARE_VOL_SENS does some bookkeeping to ensure that the volume
% conductor model and the sensor array are ready for subsequent forward
% leadfield computations. It takes care of some pre-computations that can
% be done efficiently prior to the leadfield calculations.
%
% Use as
%   [vol, sens] = ft_prepare_vol_sens(vol, sens, ...)
% with input arguments
%   sens   structure with gradiometer or electrode definition
%   vol    structure with volume conductor definition
%
% The vol structure represents a volume conductor model, its contents
% depend on the type of model. The sens structure represents a sensor
% array, i.e. EEG electrodes or MEG gradiometers.
%
% Additional options should be specified in key-value pairs and can be
%   'channel'    cell-array with strings (default = 'all')
%   'order'      number, for single shell "Nolte" model (default = 10)
%
% The detailed behaviour of this function depends on whether the input
% consists of EEG or MEG and furthermoree depends on the type of volume
% conductor model:
% - in case of EEG single and concentric sphere models, the electrodes are
%   projected onto the skin surface.
% - in case of EEG boundary element models, the electrodes are projected on
%   the surface and a blilinear interpoaltion matrix from vertices to
%   electrodes is computed.
% - in case of MEG and a localspheres model, a local sphere is determined
%   for each coil in the gradiometer definition.
%  - in case of MEG with a singleshell Nolte model, the volume conduction
%    model is initialized
% In any case channel selection and reordering will be done. The channel
% order returned by this function corresponds to the order in the 'channel'
% option, or if not specified, to the order in the input sensor array.
%
% See also FT_COMPUTE_LEADFIELD, FT_READ_VOL, FT_READ_SENS, FT_TRANSFORM_VOL,
% FT_TRANSFORM_SENS

% Copyright (C) 2004-2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_prepare_vol_sens.m 8997 2013-12-09 13:23:21Z roboos $

% get the optional input arguments
% fileformat = ft_getopt(varargin, 'fileformat');


channel = ea_ft_getopt(varargin, 'channel', sens.label);   % cell-array with channel labels, default is all
order   = ea_ft_getopt(varargin, 'order', 10);             % order of expansion for Nolte method; 10 should be enough for real applications; in simulations it makes sense to go higher

% ensure that the sensor description is up-to-date (Aug 2011)
%sens = ft_datatype_sens(sens);



% ensure that the volume conduction description is up-to-date (Jul 2012)
%vol = ft_datatype_headmodel(vol);

% % determine the skin compartment
% if ~isfield(vol, 'skin_surface')
%   if isfield(vol, 'bnd')
%     vol.skin_surface   = find_outermost_boundary(vol.bnd);
%   elseif isfield(vol, 'r') && length(vol.r)<=4
%     [dum, vol.skin_surface] = max(vol.r);
%   end
% end
%
% % determine the inner_skull_surface compartment
% if ~isfield(vol, 'inner_skull_surface')
%   if isfield(vol, 'bnd')
%     vol.inner_skull_surface  = find_innermost_boundary(vol.bnd);
%   elseif isfield(vol, 'r') && length(vol.r)<=4
%     [dum, vol.inner_skull_surface] = min(vol.r);
%   end
% end

% % otherwise the voltype assignment to an empty struct below won't work
% if isempty(vol)
%   vol = [];
% end

% this makes them easier to recognise
sens.type = ea_ft_senstype(sens);
vol.type  = ea_ft_voltype(vol);

if isfield(vol, 'unit') && isfield(sens, 'unit') && ~strcmp(vol.unit, sens.unit)
    error('inconsistency in the units of the volume conductor and the sensor array');
end

% whole meg part removed from original ft_ file.

  % the electrodes are used, the channel positions are not relevant any more
  % channel positinos need to be recomputed after projecting the electrodes on the skin
  if isfield(sens, 'chanpos'); sens = rmfield(sens, 'chanpos'); end

  % select the desired channels from the electrode array
  % order them according to the users specification
  [selchan, selsens] = match_str(channel, sens.label);
  Nchans = length(sens.label);

  sens.label     = sens.label(selsens);
  try, sens.chantype  = sens.chantype(selsens); end;
  try, sens.chanunit  = sens.chanunit(selsens); end;

  if isfield(sens, 'tra')
    % first only modify the linear combination of electrodes into channels
    sens.tra     = sens.tra(selsens,:);
    % subsequently remove the electrodes that do not contribute to any channel output
    selelec      = any(sens.tra~=0,1);
    sens.elecpos = sens.elecpos(selelec,:);
    sens.tra     = sens.tra(:,selelec);
  else
    % the electrodes and channels are identical
    sens.elecpos = sens.elecpos(selsens,:);
  end

%   switch ea_ft_voltype(vol)
%
%     case 'simbio' - always simbio in our case..

      % extract the outer surface
      bnd = ea_mesh2edge(vol);
      for j=1:length(sens.label)
        d = bsxfun(@minus, bnd.pnt, sens.elecpos(j,:));
        [d, i] = min(sum(d.^2, 2));
        % replace the position of each electrode by the closest vertex
        sens.elecpos(j,:) = bnd.pnt(i,:);
      end

      vol.transfer = ea_sb_transfer(vol,sens);
  % end

  % FIXME this needs careful thought to ensure that the average referencing which is now done here and there, and that the linear interpolation in case of BEM are all dealt with consistently
  % % always ensure that there is a linear transfer matrix for
  % % rereferencing the EEG potential
  % if ~isfield(sens, 'tra');
  %   sens.tra = eye(length(sens.label));
  % end

  % update the channel positions as the electrodes were projected to the skin surface
  [pos, ~, lab] = ea_channelposition(sens);
  [selsens, selpos] = match_str(sens.label, lab);
  sens.chanpos = nan(length(sens.label),3);
  sens.chanpos(selsens,:) = pos(selpos,:);


if isfield(sens, 'tra')
  if issparse(sens.tra) && size(sens.tra, 1)==1
    % this multiplication would result in a sparse leadfield, which is not what we want
    % the effect can be demonstrated as sparse(1)*rand(1,10), see also http://bugzilla.fcdonders.nl/show_bug.cgi?id=1169#c7
    sens.tra = full(sens.tra);
  elseif ~issparse(sens.tra) && size(sens.tra, 1)>1
    % the multiplication of the "sensor" leadfield (electrode or coil) with the tra matrix to get the "channel" leadfield
    % is faster for most cases if the pre-multiplying weighting matrix is made sparse
    sens.tra = sparse(sens.tra);
  end
end

function [transfer] = ea_sb_transfer(vol,sens)

% SB_TRANSFER
%
% INPUT:
% vol.wf.nd,vol.wf.el: nodes and elements of head grid (nodes in mm)
% vol.wf.field: vector assigning a conductivity to teach element (in S/mm)
% elc: vector containing the electrode positions (in mm)
%
% OUTPUT:
% vol.transfer: transfer matrix
%
% $Id: sb_transfer.m 8776 2013-11-14 09:04:48Z roboos $

%--------------------------------------------------------------------------
%TODO: tissuecond: vector with the conductivities of the respective labels,
%right now not used
%--------------------------------------------------------------------------
%TODO: automated setting of conductivities?
%disp('Setting conductivities...')
%vol.wf.field = sb_set_cond(vol.wf.labels,tissuecond); %gucken, wo tissuecond herkommt...
%--------------------------------------------------------------------------
%TODO: this probably needs to be expanded so that the node is a surface
%node if no ECoG-flag is set - solution: find a surface first...
%--------------------------------------------------------------------------

disp('Find electrode positions...')
vol.elecnodes = ea_sb_find_elec(vol,sens);
%calculate transfermatrix
disp('Calculate transfer matrix...')
transfer = zeros(length(vol.elecnodes),size(vol.pos,1));
%--------------------------------------------------------------------------
%TODO: add possibility for parallel computation
%matlabpool local 2;

ea_dispercent(0,'Solving forward model');
for i=2:length(vol.elecnodes)
    ea_dispercent(i/length(vol.elecnodes));
    %str = ['Electrode ',num2str(i),' of ',num2str(size(vol.elecnodes,1))];
    %disp(str)
    vecb = zeros(size(vol.stiff,1),1);
    vecb(vol.elecnodes(i)) = 1;
    transfer(i,:) = ea_sb_calc_vecx(vol.stiff,vecb,vol.elecnodes(1));
    clear vecb;
end
ea_dispercent(1,'end');

function diri = ea_sb_find_elec(vol,sens)

% SB_FIND_ELEC
%
% $Id: sb_find_elec.m 8776 2013-11-14 09:04:48Z roboos $

diri = zeros(size(sens.elecpos,1),1);
for i=1:size(sens.elecpos,1)
    [dist, diri(i)] = min(sum(bsxfun(@minus,vol.pos,sens.elecpos(i,:)).^2,2));
end

function vecx = ea_sb_calc_vecx(stiff,vecb,ref)

% SB_CALC_VECX
%
% $Id: sb_calc_vecx.m 8776 2013-11-14 09:04:48Z roboos $

vecdi = zeros(size(stiff,1),1);
vecdi(ref) = 1;
vecva = zeros(size(stiff,1),1);
[stiff, vecb] = ea_sb_set_bndcon(stiff,vecb,vecdi,vecva);
clear vecdi, vecva;
vecx = ea_sb_solve(stiff,vecb);

function [stiff,rhs] = ea_sb_set_bndcon(stiff,rhs,dirinode,dirival)

% SB_SET_BNDCON
%
% $Id: sb_set_bndcon.m 8776 2013-11-14 09:04:48Z roboos $

dia = diag(stiff);
stiff = stiff - diag(dia);
[indexi indexj s] = find(stiff);
clear stiff;
dind = find(dirinode > 0);
indi = find(ismember(indexi,dind));
indj = find(~ismember(indexj,dind));
indij = intersect(indi,indj);
rhs(indexj(indij)) = rhs(indexj(indij)) + dirival(indexi(indij)).*s(indij);
s(indi) = 0;
dia(indexi(indi)) = 1;
rhs(indexi(indi)) = -dirival(indexi(indi));
indij = find(ismember(indexj,dind)&~ismember(indexi,dind));
rhs(indexi(indij)) = rhs(indexi(indij)) + dirival(indexj(indij)).*s(indij);
s(indij) = 0;
stiff = sparse(indexi,indexj,s,length(dia),length(dia));
stiff = stiff + diag(dia);

function x = ea_sb_solve(sysmat,vecb)

% SB_SOLVE
%
% $Id: sb_solve.m 8776 2013-11-14 09:04:48Z roboos $

%scalen
%disp('Scaling stiffness matrix...')
dkond = 1./(sqrt(diag(sysmat)));
vecb = vecb.*dkond;
[indexi, indexj, s] = find(sysmat);
sys_size = size(sysmat,1);
clear sysmat;
s = (s.*dkond(indexi)).*dkond(indexj);
s(1) = 1;
%disp('Preconditioning...')
L = sparse(indexi,indexj,s,sys_size,sys_size,length(s));
%partch
try
    L = ichol(L);
catch
    disp('Could not compute incomplete Cholesky-decompositon. Rescaling stiffness matrix...')
    alpha = 0.5d-6;
    alpha = alpha*8.d0;
    alpha = 1 / (alpha + 1);
    s = alpha*s;
    dia = find(indexi == indexj);
    s(dia) = (1./alpha)*s(dia);
    s(dia) = sqrt(s(dia));
    s(1) = 1;
    L = sparse(indexi,indexj,s,sys_size,sys_size,length(s));
    clear dia;
    L = ichol(L);
end
%startvektor
%disp('Finding startvector...')
vecb_ = L \ (-vecb);
vecx = L' \ vecb_;
clear vecb_;
%sonstiges
sysmat = sparse(indexi,indexj,s,sys_size,sys_size,length(s));
sysmat = sysmat + sysmat' - sparse(1:sys_size,1:sys_size,diag(sysmat),sys_size,sys_size,sys_size);
clear indexi indexj s;
%dprod = sysmat * vecx;
%loesen
%disp('Solving equation system...')
[~,x]=evalc('pcg(sysmat,vecb,10e-9,5000,L,L'',vecx)');
%fl
%rr
%it
%rescal
x = x.*dkond;

function [newbnd] = ea_mesh2edge(bnd)

% MESH2EDGE finds the edge lines from a triangulated mesh or the edge surfaces
% from a tetrahedral or hexahedral mesh.
%
% Use as
%   [bnd] = mesh2edge(bnd)

% Copyright (C) 2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: mesh2edge.m 8776 2013-11-14 09:04:48Z roboos $

if isfield(bnd, 'tri')
  % make a list of all edges
  edge1 = bnd.tri(:, [1 2]);
  edge2 = bnd.tri(:, [2 3]);
  edge3 = bnd.tri(:, [3 1]);
  edge = cat(1, edge1, edge2, edge3);

elseif isfield(bnd, 'tet')
  % make a list of all triangles that form the tetraheder
  tri1 = bnd.tet(:, [1 2 3]);
  tri2 = bnd.tet(:, [2 3 4]);
  tri3 = bnd.tet(:, [3 4 1]);
  tri4 = bnd.tet(:, [4 1 2]);
  edge = cat(1, tri1, tri2, tri3, tri4);

elseif isfield(bnd, 'hex')
  % make a list of all "squares" that form the cube/hexaheder
  % FIXME should be checked, this is impossible without a drawing
  square1 = bnd.hex(:, [1 2 3 4]);
  square2 = bnd.hex(:, [5 6 7 8]);
  square3 = bnd.hex(:, [1 2 6 5]);
  square4 = bnd.hex(:, [2 3 7 6]);
  square5 = bnd.hex(:, [3 4 8 7]);
  square6 = bnd.hex(:, [4 1 5 8]);
  edge = cat(1, square1, square2, square3, square4, square5, square6);

end % isfield(bnd)

% soort all polygons in the same direction
% keep the original as "edge" and the sorted one as "sedge"
sedge = sort(edge, 2);

% % find the edges that are not shared -> count the number of occurences
% n = size(sedge,1);
% occurences = ones(n,1);
% for i=1:n
%   for j=(i+1):n
%     if all(sedge(i,:)==sedge(j,:))
%       occurences(i) = occurences(i)+1;
%       occurences(j) = occurences(j)+1;
%     end
%   end
% end
%
% % make the selection in the original, not the sorted version of the edges
% % otherwise the orientation of the edges might get flipped
% edge = edge(occurences==1,:);

% find the edges that are not shared
indx = ea_findsingleoccurringrows(sedge);
edge = edge(indx, :);

if ~isfield(bnd, 'pnt') && isfield(bnd, 'pos')
  bnd.pnt = bnd.pos;
end

% the naming of the output edges depends on what they represent
newbnd.pnt  = bnd.pnt;
if isfield(bnd, 'tri')
  newbnd.line = edge;
elseif isfield(bnd, 'tet')
  newbnd.tri = edge;
elseif isfield(bnd, 'hex')
  newbnd.poly = edge;
end

function indx = ea_findsingleoccurringrows(X)
[X, indx] = sortrows(X);
sel  = any(diff([X(1,:)-1; X],1),2) & any(diff([X; X(end,:)+1],1),2);
indx = indx(sel);

function [type] = ea_ft_voltype(vol, desired)

% FT_VOLTYPE determines the type of volume conduction model of the head
%
% Use as
%   [type] = ft_voltype(vol)
% to get a string describing the type, or
%   [flag] = ft_voltype(vol, desired)
% to get a boolean value.
%
% For EEG the following volume conduction models are recognized
%   singlesphere       analytical single sphere model
%   concentricspheres  analytical concentric sphere model with up to 4 spheres
%   halfspace          infinite homogenous medium on one side, vacuum on the other
%   openmeeg           boundary element method, based on the OpenMEEG software
%   bemcp              boundary element method, based on the implementation from Christophe Phillips
%   dipoli             boundary element method, based on the implementation from Thom Oostendorp
%   asa                boundary element method, based on the (commercial) ASA software
%   simbio             finite element method, based on the SimBio software
%   fns                finite difference method, based on the FNS software
%   interpolate        interpolate the potential based on pre-computed leadfields
%
% and for MEG the following volume conduction models are recognized
%   singlesphere       analytical single sphere model
%   localspheres       local spheres model for MEG, one sphere per channel
%   singleshell        realisically shaped single shell approximation, based on the implementation from Guido Nolte
%   infinite           magnetic dipole in an infinite vacuum
%   interpolate        interpolate the potential based on pre-computed leadfields
%
% See also FT_COMPUTE_LEADFIELD, FT_READ_VOL, FT_HEADMODEL_BEMCP,
% FT_HEADMODEL_ASA, FT_HEADMODEL_DIPOLI, FT_HEADMODEL_SIMBIO,
% FT_HEADMODEL_FNS, FT_HEADMODEL_HALFSPACE, FT_HEADMODEL_INFINITE,
% FT_HEADMODEL_OPENMEEG, FT_HEADMODEL_SINGLESPHERE,
% FT_HEADMODEL_CONCENTRICSPHERES, FT_HEADMODEL_LOCALSPHERES,
% FT_HEADMODEL_SINGLESHELL, FT_HEADMODEL_INTERPOLATE

% Copyright (C) 2007-2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_voltype.m 10012 2014-12-03 09:14:18Z roboos $

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

if iscell(vol) && numel(vol)<4
  % this might represent combined EEG, ECoG and/or MEG
  type = cell(size(vol));
  if nargin<2
    desired = cell(size(vol)); % empty elements
  end
  for i=1:numel(vol)
    type{i} = ea_ft_voltype(vol{i}, desired{i});
  end
  return
end

if nargin<2
  % ensure that all input arguments are defined
  desired = [];
end

current_argin = {vol, desired};
if isequal(current_argin, previous_argin)
  % don't do the type detection again, but return the previous values from cache
  type = previous_argout{1};
  return
end

if isfield(vol, 'type') && ~(ft_datatype(vol, 'grad') || ft_datatype(vol, 'sens')) % grad and sens also contain .type fields
  % preferably the structure specifies its own type
  type = vol.type;

elseif isfield(vol, 'r') && numel(vol.r)==1 && ~isfield(vol, 'label')
  type = 'singlesphere';

elseif isfield(vol, 'r') && isfield(vol, 'o') && isfield(vol, 'label')
  % this is before the spheres have been assigned to the coils
  % and every sphere is still associated with a channel
  type = 'localspheres';

elseif isfield(vol, 'r') && isfield(vol, 'o') && size(vol.r,1)==size(vol.o,1) && size(vol.r,1)>4
  % this is after the spheres have been assigned to the coils
  % note that this one is easy to confuse with the concentric one
  type = 'localspheres';

elseif isfield(vol, 'r') && numel(vol.r)>=2 && ~isfield(vol, 'label')
  type = 'concentricspheres';

elseif isfield(vol, 'bnd') && isfield(vol, 'mat')
  type = 'bem'; % it could be dipoli, asa, bemcp or openmeeg

elseif isfield(vol, 'bnd') && isfield(vol, 'forwpar')
  type = 'singleshell';

elseif isfield(vol, 'bnd') && numel(vol.bnd)==1
  type = 'singleshell';

elseif isempty(vol) || (isstruct(vol) && isequal(fieldnames(vol), {'unit'}))
  % it is empty, or only contains a specification of geometrical units
  type = 'infinite';

else
  type = 'unknown';

end % if isfield(vol, 'type')

if ~isempty(desired)
  % return a boolean flag
  switch desired
    case 'bem'
      type = any(strcmp(type, {'bem', 'dipoli', 'asa', 'bemcp', 'openmeeg'}));
    otherwise
      type = any(strcmp(type, desired));
  end % switch desired
end % determine the correspondence to the desired type

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout  = {type};
previous_argin  = current_argin;
previous_argout = current_argout;

return % voltype main()

function [bnd, cfg] = ea_ft_prepare_mesh(cfg, mri)

% FT_PREPARE_MESH creates a triangulated surface mesh for the volume
% conduction model. The mesh can either be selected manually from raw
% mri data or can be generated starting from a segmented volume
% information stored in the mri structure. The result is a bnd
% structure which contains the information about all segmented surfaces
% related to mri and are expressed in world coordinates.
%
% Use as
%   bnd = ft_prepare_mesh(cfg, mri)
%   bnd = ft_prepare_mesh(cfg, seg)
%
% Configuration options:
%   cfg.method      = string, can be 'interactive', 'projectmesh', 'iso2mesh', 'isosurface',
%                     'headshape', 'hexahedral', 'tetrahedral'
%   cfg.tissue      = cell-array with tissue types or numeric vector with integer values
%   cfg.numvertices = numeric vector, should have same number of elements as cfg.tissue
%   cfg.downsample  = integer number (default = 1, i.e. no downsampling), see FT_VOLUMEDOWNSAMPLE
%   cfg.headshape   = (optional) a filename containing headshape, a Nx3 matrix with surface
%                     points, or a structure with a single or multiple boundaries
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% Example
%   mri             = ft_read_mri('Subject01.mri');
%
%   cfg             = [];
%   cfg.output      = {'scalp', 'skull', 'brain'};
%   segmentation    = ft_volumesegment(cfg, mri);
%
%   cfg             = [];
%   cfg.tissue      = {'scalp', 'skull', 'brain'};
%   cfg.numvertices = [800, 1600, 2400];
%   bnd             = ft_prepare_mesh(cfg, segmentation);
%
% See also FT_VOLUMESEGMENT, FT_PREPARE_HEADMODEL, FT_PLOT_MESH

% Undocumented functionality: at this moment it allows for either
%   bnd = ft_prepare_mesh(cfg)             or
%   bnd = ft_prepare_mesh(cfg, headmodel)
% but more consistent would be to specify a volume conduction model with
%   cfg.vol           = structure with volume conduction model, see FT_PREPARE_HEADMODEL
%   cfg.headshape     = name of file containing the volume conduction model, see FT_READ_VOL
%
% Undocumented options, I have no clue why they exist
%   cfg.method = {'singlesphere' 'concentricspheres' 'localspheres'}

% Copyrights (C) 2009-2012, Robert Oostenveld & Cristiano Micheli
% Copyrights (C) 2012-2013, Robert Oostenveld & Lilla Magyari
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_prepare_mesh.m 9654 2014-06-21 07:24:04Z roboos $

revision = '$Id: ft_prepare_mesh.m 9654 2014-06-21 07:24:04Z roboos $';

% do the general setup of the function
ea_ft_defaults
ea_ft_preamble init
ea_ft_preamble provenance
ea_ft_preamble trackconfig
ea_ft_preamble debug
ea_ft_preamble loadvar mri

% the abort variable is set to true or false in ft_preamble_init
% if abort
%   return
% end

% we cannot use nargin, because the data might have been loaded from cfg.inputfile
hasdata = exist('mri', 'var');

% check if the input cfg is valid for this function
%cfg = ft_checkconfig(cfg, 'forbidden', {'numcompartments', 'outputfile', 'sourceunits', 'mriunits'});

% get the options
cfg.downsample  = ea_ft_getopt(cfg, 'downsample', 1); % default is no downsampling
cfg.numvertices = ea_ft_getopt(cfg, 'numvertices');   % no default

% This was changed on 3 December 2013, this backward compatibility can be removed in 6 months from now.
if isfield(cfg, 'interactive')
  if strcmp(cfg.interactive, 'yes')
    warning('please specify cfg.method=''interactive'' instead of cfg.interactive=''yes''');
    cfg.method = 'interactive';
  end
  cfg = rmfield(cfg, 'interactive');
end

% This was changed on 3 December 2013, it makes sense to keep it like this on the
% long term (previously there was no explicit use of cfg.method, now there is).
% Translate the input options in the appropriate cfg.method.
if ~isfield(cfg, 'method')
  if isfield(cfg, 'headshape') && ~isempty(cfg.headshape)
    warning('please specify cfg.method=''headshape''');
    cfg.method = 'headshape';
  elseif hasdata && ~strcmp(ea_ft_voltype(mri), 'unknown')
    % the input is a spherical volume conduction model
    cfg.method = ea_ft_voltype(mri);
  elseif hasdata
    warning('please specify cfg.method=''projectmesh'', ''iso2mesh'' or ''isosurface''');
    warning('using ''projectmesh'' as default');
    cfg.method = 'projectmesh';
  end
end

%% the following does not apply in our case
% if hasdata && cfg.downsample~=1
%   % optionally downsample the anatomical volume and/or tissue segmentations
%   tmpcfg = keepfields(cfg, {'downsample'});
%   mri = ft_volumedownsample(tmpcfg, mri);
%   [cfg, mri] = rollback_provenance(cfg, mri);
% end

switch cfg.method
  % deleted other methods from original fieldtrip function..
  case 'hexahedral'
    % the MRI is assumed to contain a segmentation
    % call the corresponding helper function
    bnd = ea_prepare_mesh_hexahedral(cfg, mri);

end

% copy the geometrical units from the input to the output
if ~isfield(bnd, 'unit') && hasdata && isfield(mri, 'unit')
  for i=1:numel(bnd)
    bnd(i).unit = mri.unit;
  end
elseif ~isfield(bnd, 'unit')
  bnd = ft_convert_units(bnd);
end

% do the general cleanup and bookkeeping at the end of the function
ea_ft_postamble debug
ea_ft_postamble trackconfig
ea_ft_postamble provenance
ea_ft_postamble previous mri
ea_ft_postamble history bnd

function [pnt, ori, lab] = ea_channelposition(sens)

% CHANNELPOSITION computes the channel positions and orientations from the
% coils or electrodes
%
% Use as
%   [pos, ori, lab] = channelposition(sens)
% where sens is an electrode or gradiometer array.
%
% See also FT_DATATYPE_SENS

% Copyright (C) 2009-2012, Robert Oostenveld & Vladimir Litvak
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: channelposition.m 8920 2013-11-29 13:20:51Z roboos $

% remove the balancing from the sensor definition, e.g. planar gradients, 3rd-order gradients, PCA-cleaned data or ICA projections
if isfield(sens, 'balance') && ~strcmp(sens.balance.current, 'none')
  sens = ea_undobalancing(sens);
end

% keep it backward compatible with sensor definitions prior to 2011v1 (see ft_datatype_sens), which have pnt/ori instead of coilpos/coilori.
if isfield(sens, 'ori')
  sens.coilori = sens.ori;
  sens = rmfield(sens, 'ori');
end
if isfield(sens, 'pnt')
  sens.coilpos = sens.pnt;
  sens = rmfield(sens, 'pnt');
end

% treat all sensor arrays similar, i.e. as gradiometer systems
if     ~isfield(sens, 'coilori') && isfield(sens, 'coilpos')
  sens.coilori = nan(size(sens.coilpos));
elseif ~isfield(sens, 'coilori') && isfield(sens, 'elecpos')
  sens.coilori = nan(size(sens.elecpos));
end

switch ft_senstype(sens)
  case {'ctf64', 'ctf151', 'ctf275' 'bti148', 'bti248', 'bti248grad', 'itab28', 'itab153', 'yokogawa64', 'yokogawa160', 'babysquid74'}
    % the following code applies to systems with only axial gradiometers or magnetometers

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % do the MEG sensors first
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    sensorig   = sens;
    sel        = ft_chantype(sens, 'meg');
    sens.label = sens.label(sel);
    sens.tra   = sens.tra(sel,:);

    % subsequently remove the unused coils
    used = any(abs(sens.tra)>0.0001, 1);  % allow a little bit of rounding-off error
    sens.coilpos = sens.coilpos(used,:);
    sens.coilori = sens.coilori(used,:);
    sens.tra     = sens.tra(:,used);

    % compute distances from the center of the helmet
    center = mean(sens.coilpos);
    dist   = sqrt(sum((sens.coilpos - repmat(center, size(sens.coilpos, 1), 1)).^2, 2));

    % put the corresponding distances instead of non-zero tra entries
    maxval = repmat(max(abs(sens.tra),[],2), [1 size(sens.tra,2)]);
    maxval = min(maxval, ones(size(maxval))); %a value > 1 sometimes leads to problems; this is an empirical fix
    dist = (abs(sens.tra)>0.7.*maxval).*repmat(dist', size(sens.tra, 1), 1);

    % for the occasional case where there are nans: -> 0's will be
    % converted to inf anyhow
    dist(isnan(dist)) = 0;

    % put infs instead of the zero entries
    dist(~dist) = inf;

    % use the matrix to find coils with minimal distance to the center,
    % i.e. the bottom coil in the case of axial gradiometers
    % this only works for a full-rank unbalanced tra-matrix

    numcoils = sum(isfinite(dist),2);

    if all(numcoils==numcoils(1))
      % add the additional constraint that coils cannot be used twice,
      % i.e. for the position of 2 channels. A row of the dist matrix can end
      % up with more than 1 (magnetometer array) or 2 (axial gradiometer array)
      % non-zero entries when the input grad structure is rank-reduced
      % FIXME: I don't know whether this works for a vector-gradiometer
      % system. It also does not work when the system has mixed gradiometers
      % and magnetometers

      % use the magic that Jan-Mathijs implemented
      tmp      = mode(numcoils);
      niter    = 0;
      while ~all(numcoils==tmp)
        niter    = niter + 1;
        selmode  = find(numcoils==tmp);
        selrest  = setdiff((1:size(dist,1))', selmode);
        dist(selrest,sum(~isinf(dist(selmode,:)))>0) = inf;
        numcoils = sum(isfinite(dist),2);
        if niter>500
          error('Failed to extract the positions of the channels. This is most likely due to the balancing matrix being rank deficient. Please replace data.grad with the original grad-structure obtained after reading the header.');
        end
      end
    else
      % assume that the solution is not so hard and just determine the bottom coil
    end

    [junk, ind] = min(dist, [], 2);

    lab(sel)   = sens.label;
    pnt(sel,:) = sens.coilpos(ind, :);
    ori(sel,:) = sens.coilori(ind, :);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % then do the references
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sens = sensorig;
    sel  = ft_chantype(sens, 'ref');

    sens.label = sens.label(sel);
    sens.tra   = sens.tra(sel,:);

    % subsequently remove the unused coils
    used = any(abs(sens.tra)>0.0001, 1);  % allow a little bit of rounding-off error
    sens.coilpos = sens.coilpos(used,:);
    sens.coilori = sens.coilori(used,:);
    sens.tra = sens.tra(:,used);

    [nchan, ncoil] = size(sens.tra);
    refpnt = zeros(nchan,3);
    refori = zeros(nchan,3); % FIXME not sure whether this will work
    for i=1:nchan
      weight = abs(sens.tra(i,:));
      weight = weight ./ sum(weight);
      refpnt(i,:) = weight * sens.coilpos;
      refori(i,:) = weight * sens.coilori;
    end
    reflab = sens.label;

    lab(sel)   = reflab;
    pnt(sel,:) = refpnt;
    ori(sel,:) = refori;

    sens = sensorig;

  case {'ctf64_planar', 'ctf151_planar', 'ctf275_planar', 'bti148_planar', 'bti248_planar', 'bti248grad_planar', 'itab28_planar', 'itab153_planar', 'yokogawa64_planar', 'yokogawa160_planar'}
    % create a list with planar channel names
    chan = {};
    for i=1:length(sens.label)
      if ~isempty(findstr(sens.label{i}, '_dH')) || ~isempty(findstr(sens.label{i}, '_dV'))
        chan{i} = sens.label{i}(1:(end-3));
      end
    end
    chan = unique(chan);
    % find the matching channel-duplets
    ind = [];
    lab = {};
    for i=1:length(chan)
      ch1 =  [chan{i} '_dH'];
      ch2 =  [chan{i} '_dV'];
      sel = match_str(sens.label, {ch1, ch2});
      if length(sel)==2
        ind = [ind; i];
        lab(i,:) = {ch1, ch2};
        meanpnt1 = mean(sens.coilpos(abs(sens.tra(sel(1),:))>0.5, :), 1);
        meanpnt2 = mean(sens.coilpos(abs(sens.tra(sel(2),:))>0.5, :), 1);
        pnt(i,:) = mean([meanpnt1; meanpnt2], 1);
        meanori1 = mean(sens.coilori(abs(sens.tra(sel(1),:))>0.5, :), 1);
        meanori2 = mean(sens.coilori(abs(sens.tra(sel(2),:))>0.5, :), 1);
        ori(i,:) = mean([meanori1; meanori2], 1);
      end
    end
    lab = lab(ind,:);
    pnt = pnt(ind,:);
    ori = ori(ind,:);

  case 'neuromag122'
    % find the matching channel-duplets
    ind = [];
    lab = {};
    for i=1:2:140
      % first try MEG channel labels with a space
      ch1 = sprintf('MEG %03d', i);
      ch2 = sprintf('MEG %03d', i+1);
      sel = match_str(sens.label, {ch1, ch2});
      % then try MEG channel labels without a space
      if (length(sel)~=2)
        ch1 = sprintf('MEG%03d', i);
        ch2 = sprintf('MEG%03d', i+1);
        sel = match_str(sens.label, {ch1, ch2});
      end
      % then try to determine the channel locations
      if (length(sel)==2)
        ind = [ind; i];
        lab(i,:) = {ch1, ch2};
        meanpnt1 = mean(sens.coilpos(abs(sens.tra(sel(1),:))>0.5,:), 1);
        meanpnt2 = mean(sens.coilpos(abs(sens.tra(sel(2),:))>0.5,:), 1);
        pnt(i,:) = mean([meanpnt1; meanpnt2], 1);
        meanori1 = mean(sens.coilori(abs(sens.tra(sel(1),:))>0.5,:), 1);
        meanori2 = mean(sens.coilori(abs(sens.tra(sel(2),:))>0.5,:), 1);
        ori(i,:) = mean([meanori1; meanori2], 1);
      end
    end
    lab = lab(ind,:);
    pnt = pnt(ind,:);
    ori = ori(ind,:);

  case 'neuromag306'
    % find the matching channel-triplets
    ind = [];
    lab = {};
    for i=1:300
      % first try MEG channel labels with a space
      ch1 = sprintf('MEG %03d1', i);
      ch2 = sprintf('MEG %03d2', i);
      ch3 = sprintf('MEG %03d3', i);
      [sel1, sel2] = match_str(sens.label, {ch1, ch2, ch3});
      % the try MEG channels without a space
      if isempty(sel1)
        ch1 = sprintf('MEG%03d1', i);
        ch2 = sprintf('MEG%03d2', i);
        ch3 = sprintf('MEG%03d3', i);
        [sel1, sel2] = match_str(sens.label, {ch1, ch2, ch3});
      end
      % then try to determine the channel locations
      if (~isempty(sel1) && length(sel1)<=3)
        ind = [ind; i];
        lab(i,sel2) = sens.label(sel1)';
        meanpnt  = [];
        meanori  = [];
        for j = 1:length(sel1)
          meanpnt  = [meanpnt; mean(sens.coilpos(abs(sens.tra(sel1(j),:))>0.5,:), 1)];
          meanori  = [meanori; mean(sens.coilori(abs(sens.tra(sel1(j),:))>0.5,:), 1)];
        end
        pnt(i,:) = mean(meanpnt, 1);
        ori(i,:) = mean(meanori, 1);
      end
    end
    lab = lab(ind,:);
    pnt = pnt(ind,:);
    ori = ori(ind,:);

  otherwise
    % compute the position for each gradiometer or electrode
    nchan = length(sens.label);
    if isfield(sens, 'elecpos')
      nelec = size(sens.elecpos,1); % these are the electrodes
    elseif isfield(sens, 'coilpos')
      ncoil = size(sens.coilpos,1); % these are the coils
    end

    if ~isfield(sens, 'tra') && isfield(sens, 'elecpos') && nchan==nelec
      % there is one electrode per channel, which means that the channel position is identical to the electrode position
      pnt = sens.elecpos;
      ori = nan(size(pnt));
      lab = sens.label;

    elseif isfield(sens, 'tra') && isfield(sens, 'elecpos') && isequal(sens.tra, eye(nelec))
      % there is one electrode per channel, which means that the channel position is identical to the electrode position
      pnt = sens.elecpos;
      ori = nan(size(pnt));
      lab = sens.label;

    elseif isfield(sens, 'tra') && isfield(sens, 'elecpos') && isequal(sens.tra, eye(nelec)-1/nelec)
      % there is one electrode per channel, channels are average referenced
      pnt = sens.elecpos;
      ori = nan(size(pnt));
      lab = sens.label;

    elseif ~isfield(sens, 'tra') && isfield(sens, 'coilpos') && nchan==ncoil
      % there is one coil per channel, which means that the channel position is identical to the coil position
      pnt = sens.coilpos;
      ori = sens.coilori;
      lab = sens.label;

    elseif isfield(sens, 'tra')
      % each channel depends on multiple sensors (electrodes or coils), compute a weighted position for the channel
      % for MEG gradiometer channels this means that the position is in between the two coils
      % for bipolar EEG channels this means that the position is in between the two electrodes
      pnt = nan(nchan,3);
      ori = nan(nchan,3);
      if isfield(sens, 'coilpos')
        for i=1:nchan
          weight = abs(sens.tra(i,:));
          weight = weight ./ sum(weight);
          pnt(i,:) = weight * sens.coilpos;
          ori(i,:) = weight * sens.coilori;
        end
      elseif isfield(sens, 'elecpos')
        for i=1:nchan
          weight = abs(sens.tra(i,:));
          weight = weight ./ sum(weight);
          pnt(i,:) = weight * sens.elecpos;
        end
      end
      lab = sens.label;

    end

end % switch senstype

n   = size(lab,2);
% this is to fix the planar layouts, which cannot be plotted anyway
if n>1 && size(lab, 1)>1 % this is to prevent confusion when lab happens to be a row array
  pnt = repmat(pnt, n, 1);
  ori = repmat(ori, n, 1);
end

% ensure that the channel order is the same as in sens
[sel1, sel2] = match_str(sens.label, lab);
lab = lab(sel2);
pnt = pnt(sel2, :);
ori = ori(sel2, :);

% ensure that it is a row vector
lab = lab(:);

% do a sanity check on the number of positions
nchan = numel(sens.label);
if length(lab)~=nchan || size(pnt,1)~=nchan || size(ori,1)~=nchan
  warning('the positions were not determined for all channels');
end

function sens = ea_undobalancing(sens)

% UNDOBALANCING removes all balancing coefficients from the gradiometer sensor array
%
% This is used in CHANNELPOSITION, FT_PREPARE_LAYOUT, FT_SENSTYPE

while isfield(sens, 'balance') && isfield(sens.balance, 'current') && ~strcmp(sens.balance.current, 'none')
  fnames = setdiff(fieldnames(sens.balance), 'current');
  indx   = find(ismember(fnames, sens.balance.current));

  if length(indx)==1,
    % undo the synthetic gradient balancing
    fprintf('undoing the %s balancing for the gradiometer definition\n', sens.balance.current);

    % if componentanalysis was followed by rejectcomponent, the balancing matrix is rank deficient
    % leading to problems in the correct allocation of the coils to the channels
    if strcmp(sens.balance.current, 'invcomp') && strcmp(sens.balance.previous{1}, 'comp')
      tra1 = full(sens.balance.invcomp.tra);
      tra2 = full(sens.balance.comp.tra);
      tra3 = tra1;
      tmp  = tra1*tra2;
      tmp  = null(tmp); % nullspace after componentanalysis and rejectcomponent
      tmp  = tmp*tmp';  % this is the part which was removed at some point
      [ix,iy]     = match_str(sens.balance.comp.labelorg, sens.balance.invcomp.labelnew);
      tra3(iy,iy) = (eye(numel(ix))+tmp(ix,ix))*tra1(iy,iy);
      sens.balance.invcomp.tra = tra3;
      % FIXME check whether this is robust
    end

    if strcmp(sens.balance.current, 'planar')
      if isfield(sens, 'type') && contains(sens.type, '_planar')
        % remove the planar postfox from the sensor type
        sens.type = sens.type(1:(end-7));
      end
    end

    sens = ft_apply_montage(sens, sens.balance.(sens.balance.current), 'inverse', 'yes', 'keepunused', 'yes', 'warning', 'no');

    if ~isfield(sens, 'chanpos') || any(isnan(sens.chanpos(:))) || any(isnan(sens.chanori(:)))
      % this happens if the data has been component-analyzed
      % try to reconstruct the channel position and orientation
      [pos, ori, lab] = channelposition(sens);
      [sel1, sel2] = match_str(sens.label, lab);
      sens.chanpos(sel1,:) = pos(sel2,:);
      sens.chanori(sel1,:) = ori(sel2,:);
    end

  else
    warning('cannot undo %s balancing in the gradiometer definition\n', sens.balance.current);
    break
  end
end

function mesh=ea_prepare_mesh_hexahedral(cfg,mri)

% PREPARE_MESH_HEXAHEDRAL
%
% See also PREPARE_MESH_SEGMENTATION, PREPARE_MESH_MANUAL, PREPARE_MESH_HEADSHAPE
%
% Configuration options for generating a regular 3-D grid
%   cfg.tissue = cell with the names of the compartments that should be
%   meshed
%   cfg.resolution = desired resolution of the mesh (standard = 1)

% Copyrights (C) 2012-2013, Johannes Vorwerk
%
% $Id: prepare_mesh_hexahedral.m 10264 2015-03-02 10:14:28Z eelspa $

% ensure that the input is consistent with what this function expects
mri = ft_checkdata(mri, 'datatype', {'volume', 'segmentation'}, 'hasunit', 'yes');

% get the default options
cfg.tissue      = ft_getopt(cfg, 'tissue');
cfg.resolution  = ft_getopt(cfg, 'resolution');
cfg.shift       = ft_getopt(cfg, 'shift');
cfg.background  = ft_getopt(cfg, 'background');

if isempty(cfg.tissue)
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
  fn = fieldnames(mri);
    for i=1:numel(fn),if (numel(mri.(fn{i}))==prod(mri.dim))&(~strcmp(fn{i},'inside')), segfield=fn{i};end;end
  cfg.tissue=setdiff(unique(mri.(segfield)(:)),0);
end

if ischar(cfg.tissue)
  % it should either be something like {'brain', 'skull', 'scalp'}, or something like [1 2 3]
  cfg.tissue = {cfg.tissue};
end

if iscell(cfg.tissue)
  % the code below assumes that it is a probabilistic representation
  if any(strcmp(cfg.tissue, 'brain'))
    mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic', 'hasbrain', 'yes');
  else
    mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic');
  end
else
  % the code below assumes that it is an indexed representation
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
end

if isempty(cfg.resolution)
    warning('Using standard resolution 1 mm')
    cfg.resolution = 1;
end

if isempty(cfg.shift)
    warning('No node-shift selected')
    cfg.shift = 0;
elseif cfg.shift > 0.3
    warning('Node-shift should not be larger than 0.3')
    cfg.shift = 0.3;
end

if isempty(cfg.background)
    cfg.background = 0;
end

% do the mesh extraction
% this has to be adjusted for FEM!!!
if iscell(cfg.tissue)
  % this assumes that it is a probabilistic representation
  % for example {'brain', 'skull', scalp'}
  try
    temp = zeros(size(mri.(cfg.tissue{1})(:)));
    for i=1:numel(cfg.tissue)
      temp = [temp,mri.(cfg.tissue{i})(:)];
    end
    [val,seg] = max(temp,[],2);
    seg = seg - 1;
    seg = reshape(seg,mri.dim);
  catch
    error('Please specify cfg.tissue to correspond to tissue types in the segmented MRI')
  end
  tissue = cfg.tissue;
else
  % this assumes that it is an indexed representation
  % for example [3 2 1]
  seg = zeros(mri.dim);
  tissue = {};
  for i=1:numel(cfg.tissue)
    seg = seg + i*(mri.seg==cfg.tissue(i));
    if isfield(mri, 'seglabel')
      try
        tissue{i} = mri.seglabel{cfg.tissue(i)};
      catch
        error('Please specify cfg.tissue to correspond to (the name or number of) tissue types in the segmented MRI')
      end
    else
      tissue{i} = sprintf('tissue %d', i);
    end
  end
end

% reslice to desired resolution

if (cfg.resolution ~= 1)
    % this should be done like this: split seg into probabilistic, reslice
    % single compartments, take maximum values
    seg_array = [];

    seg_indices = unique(seg);

    for i=1:(length(unique(seg)))
        seg_reslice.anatomy = double(seg == (i-1));
        seg_reslice.dim = mri.dim;
        seg_reslice.transform = eye(4);
        seg_reslice.transform(1:3,4) = -ceil(mri.dim/2);

        cfg_reslice = [];
        cfg_reslice.resolution = cfg.resolution;
        cfg_reslice.dim = ceil(mri.dim/cfg.resolution);

        seg_build = ft_volumereslice(cfg_reslice,seg_reslice);

        seg_array = [seg_array,seg_build.anatomy(:)];

        clear seg_reslice;
    end

    [max_seg seg_build.seg] = max(seg_array,[],2);

    clear max_seg seg_array;

    seg_build.seg = reshape(seg_build.seg,seg_build.dim);
    seg_build.seg = seg_indices(seg_build.seg);
    seg_build.transform = mri.transform;

    clear seg_build.anatomy;
else
    seg_build.seg = seg;
    seg_build.dim = mri.dim;

    clear seg;
end

% ensure that the segmentation is binary and that there is a single contiguous region
% FIXME is this still needed when it is already binary?
%seg = volumethreshold(seg, 0.5, tissue);

ft_hastoolbox('simbio', 1);

% build the mesh

mesh = build_mesh_hexahedral(cfg,seg_build);

% converting position of meshpoints to the head coordinate system

if (cfg.resolution ~= 1)
    mesh.pnt = cfg.resolution * mesh.pnt;
end

mesh.pnt = ft_warp_apply(mri.transform,mesh.pnt,'homogeneous');

labels = mesh.labels;

clear mesh.labels;

mesh.tissue = zeros(size(labels));
numlabels = size(unique(labels),1);
mesh.tissuelabel = {};
ulabel = sort(unique(labels));
for i = 1:numlabels
  mesh.tissue(labels == ulabel(i)) = i;
  mesh.tissuelabel{i} = tissue{i};
end

function mesh = build_mesh_hexahedral(cfg,mri)

background = cfg.background;
shift = cfg.shift;
% extract number of voxels in each direction
% x_dim = mri.dim(1);
% y_dim = mri.dim(2);
% z_dim = mri.dim(3);
%labels = mri.seg;
fprintf('Dimensions of the segmentation before restriction to bounding-box: %i %i %i\n', mri.dim(1), mri.dim(2), mri.dim(3));

[bb_x, bb_y, bb_z] = ind2sub(size(mri.seg),find(mri.seg));
shift_coord = [min(bb_x) - 2,min(bb_y) - 2,min(bb_z) - 2];
bb_x = [min(bb_x), max(bb_x)];
bb_y = [min(bb_y), max(bb_y)];
bb_z = [min(bb_z), max(bb_z)];
x_dim = size(bb_x(1)-1:bb_x(2)+1,2);
y_dim = size(bb_y(1)-1:bb_y(2)+1,2);
z_dim = size(bb_z(1)-1:bb_z(2)+1,2);
labels = zeros(x_dim,y_dim,z_dim);
labels(2:(x_dim-1),2:(y_dim-1),2:(z_dim-1)) = mri.seg(bb_x(1):bb_x(2),bb_y(1):bb_y(2),bb_z(1):bb_z(2));

fprintf('Dimensions of the segmentation after restriction to bounding-box: %i %i %i\n', x_dim, y_dim, z_dim);

% create elements

mesh.hex = create_elements(x_dim,y_dim,z_dim);
fprintf('Created elements...\n' )


% create nodes

mesh.pnt = create_nodes(x_dim,y_dim,z_dim);
fprintf('Created nodes...\n' )


if(shift < 0 | shift > 0.3)
    error('Please choose a shift parameter between 0 and 0.3!');
elseif(shift > 0)

    mesh.pnt = shift_nodes(mesh.pnt,mesh.hex,labels, shift,x_dim,y_dim,z_dim);

end

%background = 1;
% delete background voxels(if desired)
if(background == 0)
mesh.hex = mesh.hex(labels ~= 0,:);
mesh.labels = labels(labels ~= 0);
else
mesh.labels = labels(:);
end


% delete unused nodes
[C, ia, ic] = unique(mesh.hex(:));
mesh.pnt = mesh.pnt(C,:,:,:);
mesh.pnt = mesh.pnt + repmat(shift_coord,size(mesh.pnt,1),1);
mesh.hex(:) = ic;

% function creating elements from a MRI-Image with the dimensions x_dim,
% y_dim, z_dim. Each voxel of the MRI-Image corresponds to one element in
% the hexahedral mesh. The numbering of the elements is as follows:
% the first x-y-plane(z == 1) is numbered by incrementing in x-direction
% first until x_dim is reached. This is done for each row in y-direction
% until y_dim is reached. Using the resulting x_dim-by-y_dim numbering as
% an overlay and adding (i-1)*(x_dim*y_dim) to the overlay while i is
% iterating through the z-dimension(until i==z_dim) we obtain a numbering
% for all the elements.
% The node-numbering is done accordingly: the bottom left node in element i
% has number i in the node-numbering. All the remaining nodes are numbered
% in the same manner as described above for the element numbering(note the
% different dimensionalities: x_dim+1 instead of x_dim etc.).
function elements = create_elements(x_dim,y_dim,z_dim)
    elements = zeros(x_dim*y_dim*z_dim,8);
    % create an offset vector for the bottom-left nodes in each element
    b = 1:((x_dim+1)*(y_dim));
    % delete the entries where the node does not correspond to an element's
    % bottom-left node(i.e. where the x-component of the node is equal to
    % x_dim+1)
    b = b(mod(b,(x_dim+1)) ~= 0);
    % repeat offset to make it fit the number of elements
    b = repmat(b,1,(z_dim));

    % create vector accounting for the offset of the nodes in z-direction
    c = fix([0:((x_dim)*(y_dim)*(z_dim)-1)]/((x_dim)*(y_dim))) * (x_dim+1)*(y_dim+1);

    % create the elements by assigning the nodes to them. entries 1 through 4
    % describe the bottom square of the hexahedron, entries 5 through 8 the top
    % square.
    elements(:,1) = b + c;
    elements(:,2) = b + c + 1;
    elements(:,3) = b + c + (x_dim+1) + 1;
    elements(:,4) = b + c + (x_dim+1);
    elements(:,5) = b + c + (x_dim + 1)*(y_dim+1);
    elements(:,6) = b + c + (x_dim + 1)*(y_dim+1) + 1;
    elements(:,7) = b + c + (x_dim + 1)*(y_dim+1) + (x_dim+1) + 1;
    elements(:,8) = b + c + (x_dim + 1)*(y_dim+1) + (x_dim+1);
    clear b;
    clear c;

% function creating the nodes and assigning coordinates in the
% [0,x_dim]x[0,y_dim]x[0,z_dim] box. for details on the node-numbering see
% comments for create_elements.
function nodes = create_nodes(x_dim, y_dim, z_dim)
    nodes = zeros(((x_dim+1)*(y_dim+1)*(z_dim + 1)),3);
    % offset vector for node coordinates
    b=[0:((x_dim+1)*(y_dim+1)*(z_dim+1)-1)];

    % assign coordinates within the box
    nodes(:,1) = mod(b,(x_dim+1));
    nodes(:,2) = mod(fix(b/(x_dim+1)),(y_dim+1));
    nodes(:,3) = fix(b/((x_dim + 1)*(y_dim+1)));

    clear b;

% function shifting the nodes
function nodes = shift_nodes(points,hex,labels, sh,x_dim,y_dim,z_dim)
    cfg = [];
    fprintf('Applying shift %f\n', sh);
    nodes = points;

    % helper vector for indexing the nodes
    b = 1:(x_dim+1)*(y_dim+1)*(z_dim+1);

    % vector which gives the corresponding element to a node(in the sense
    % of how the node-numbering was done, see comment for create_elements).
    % note that this vector still contains values for nodes which do not have a
    % corresponding element. this will be manually omitted in the finding
    % of the surrounding elements(see below).
    offset = (b - nodes(b,2)' - (nodes(b,3))'*(y_dim+1+x_dim))';
    offset(offset <= 0) = size(labels,1)+1;
    offset(offset > size(hex,1)) = size(labels,1)+1;

    % create array containing the surrounding elements for each node
    %surrounding = zeros((x_dim+1)*(y_dim+1)*(z_dim+1),8);

    % find out the surrounding of each node(if there is any)
    % find element to which the node is bottom left front
    %surrounding((nodes(b,1) < x_dim) & (nodes(b,3) < z_dim),1) = offset((nodes(b,1) < x_dim) & (nodes(b,3) < z_dim));
    % bottom right front
    %surrounding(nodes(b,1) > 0,2) = offset(nodes(b,1) > 0,1) -1;
    % bottom left back
    %surrounding(nodes(b,2) > 0,3) = offset(nodes(b,2) > 0,1) - x_dim;
    % bottom right back
    %surrounding((nodes(b,2) > 0) & (nodes(b,1) > 0),4) = offset((nodes(b,2) > 0) & (nodes(b,1) > 0),1) - (x_dim) - 1;
    % top left front
    %surrounding(nodes(b,3) > 0,5) = offset(nodes(b,3) > 0,1) - (x_dim)*(y_dim);
    % top right front
    %surrounding(nodes(b,3) > 0,6) = offset(nodes(b,3) > 0,1) - (x_dim)*(y_dim) - 1;
    % top left back
    %surrounding(nodes(b,3) > 0,7) = offset(nodes(b,3) > 0,1) - (x_dim)*(y_dim) - x_dim;
    % top right back
    %surrounding((nodes(b,3) > 0) & (nodes(b,2) > 0),8) = offset((nodes(b,3) > 0) & (nodes(b,2) > 0),1) - (x_dim)*(y_dim) - (x_dim) -1;
    %clear offset;
    %clear b;

    % those entries in the surrounding matrix which point to non-existing
    % elements(> size(hex,1) or <= 0) we overwrite with a
    % dummy element
    %surrounding(surrounding <= 0) = size(labels,1) + 1;
    %surrounding(surrounding > size(hex,1)) = size(labels,1)+1;

    % set the label of the dummy element to be zero(background)
    labels(size(labels,1)+1) = 0;

    % matrixs holding the label of each surrounding element
    %surroundinglabels = labels(surrounding);
    surroundinglabels = zeros((x_dim+1)*(y_dim+1)*(z_dim+1),8);
    surroundinglabels((nodes(b,1) < x_dim) & (nodes(b,3) < z_dim),1) = labels(offset((nodes(b,1) < x_dim) & (nodes(b,3) < z_dim)));
    surroundinglabels(nodes(b,1) > 0,2) = labels(offset(nodes(b,1) > 0,1) -1);
    surroundinglabels(nodes(b,2) > 0,3) = labels(offset(nodes(b,2) > 0,1) - x_dim);
    offsetnow = offset((nodes(b,2) > 0) & (nodes(b,1) > 0),1) - (x_dim) - 1;
    offsetnow(offsetnow <= 0) = size(labels,1)+1;
    surroundinglabels((nodes(b,2) > 0) & (nodes(b,1) > 0),4) = labels(offsetnow);
    offsetnow = offset(nodes(b,3) > 0,1) - (x_dim)*(y_dim);
    offsetnow(offsetnow <= 0) = size(labels,1)+1;
    surroundinglabels(nodes(b,3) > 0,5) = labels(offsetnow);
    offsetnow = offset(nodes(b,3) > 0,1) - (x_dim)*(y_dim) - 1;
    offsetnow(offsetnow <= 0) = size(labels,1)+1;
    surroundinglabels(nodes(b,3) > 0,6) = labels(offsetnow);
    offsetnow = offset(nodes(b,3) > 0,1) - (x_dim)*(y_dim) - x_dim;
    offsetnow(offsetnow <= 0) = size(labels,1)+1;
    surroundinglabels(nodes(b,3) > 0,7) = labels(offsetnow);
    offsetnow = offset((nodes(b,3) > 0) & (nodes(b,2) > 0),1) - (x_dim)*(y_dim) - (x_dim) -1;
    offsetnow(offsetnow <= 0) = size(labels,1)+1;
    surroundinglabels((nodes(b,3) > 0) & (nodes(b,2) > 0),8) = labels(offsetnow);

    % matrix showing how many types of each label are around a given node
    distribution = zeros(size(nodes,1), size(unique(labels),1));
    % the following assumes that the labels are 1,2,...
    for l=1:size(unique(labels),1)
        for k=1:8
            distribution(:,l) = distribution(:,l) + (surroundinglabels(:,k) == l);
        end
    end

    % fill up the last column with the amount of background labels
    distribution(:,(size(unique(labels),1))) = 8 - sum(distribution(:,1:size(unique(labels),1))');

    % how many different labels are there around each node
    distsum = sum(distribution>0,2);

    % set zeros to Inf in order to make the finding of a minimum
    % meaningful
    distribution(distribution == 0) = Inf;

    % find out which is the minority label
    [mins,minpos] = min(distribution, [],2);
    clear distribution;

    % calculate the centroid for each element
    centroids = zeros(size(hex,1),3);
    for l=1:3
        centroids(:,l) = sum(reshape(nodes(hex(:,:),l),size(hex,1),8)')'/8;
    end

    % set a dummy centroid
    centroids(size(centroids,1) +1,:) = 0;

    tbc = zeros(size(surroundinglabels));
    % helper matrix, c(i,j,k) is one when surroundinglabels(i,j) == k
    for i=1:size(unique(labels),1)+1
        c = zeros(size(surroundinglabels,2), size(unique(labels),1)+1);
        if (i == size(unique(labels),1)+1)
           c = surroundinglabels == 0;
        else
           c = surroundinglabels == i;
        end
        tbc(ismember(minpos,i) == 1,:) = c(ismember(minpos,i) == 1,:);
    end

%     % matrix that shows which elements are to be considered as minority
%     % around a given node
%     clear surroundinglabels;
%      % helper matrix, c(i,j,k) is one when surroundinglabels(i,j) == k
%     c = zeros(size(surroundinglabels,1),size(surroundinglabels,2),size(unique(labels),1)+1);
%     for i=1:size(unique(labels),1)
%         c(:,:,i) = surroundinglabels == i;
%     end
%     c(:,:,size(unique(labels),1)+1) = surroundinglabels == 0;
%
%
%     % matrix that shows which elements are to be considered as minority
%     % around a given node
%     tbc = zeros(size(surroundinglabels));
%     clear surroundinglabels;
%     for i=1:size(unique(labels),1)+1
%     tbc(ismember(minpos,i) == 1,:) = c(ismember(minpos,i) == 1,:,i);
%     end
    clear c;

    % delete cases in which we don't have a real minimum
    tbcsum = sum(tbc,2);
    tbc(tbcsum == 8,:) = 0;
    tbc(tbcsum == 4,:) = 0;
    tbcsum((distsum>2) & (mins > 1),:) = 0;
    tbcsum(tbcsum == 8) = 0;
    tbcsum(tbcsum == 4) = 0;
    tbcsum((distsum>2) & (mins > 1)) = 0;

    %surroundingconsidered = surrounding.*tbc;
    surroundingconsidered = zeros(size(tbc));
    surroundingconsidered((nodes(b,1) < x_dim) & (nodes(b,3) < z_dim),1) = offset((nodes(b,1) < x_dim) & (nodes(b,3) < z_dim)).*tbc((nodes(b,1) < x_dim) & (nodes(b,3) < z_dim),1);
    surroundingconsidered(nodes(b,1) > 0,2) = (offset(nodes(b,1) > 0,1) -1).*(tbc(nodes(b,1) > 0,2));
    surroundingconsidered(nodes(b,2) > 0,3) = (offset(nodes(b,2) > 0,1) - x_dim).*tbc(nodes(b,2) > 0,3);
    surroundingconsidered((nodes(b,2) > 0) & (nodes(b,1) > 0),4) = (offset((nodes(b,2) > 0) & (nodes(b,1) > 0),1) - (x_dim) - 1).*tbc((nodes(b,2) > 0) & (nodes(b,1) > 0),4).*tbc((nodes(b,2) > 0) & (nodes(b,1) > 0),4);
    surroundingconsidered(nodes(b,3) > 0,5) = (offset(nodes(b,3) > 0,1) - (x_dim)*(y_dim)).*tbc(nodes(b,3) > 0,5);
    surroundingconsidered(nodes(b,3) > 0,6) = (offset(nodes(b,3) > 0,1) - (x_dim)*(y_dim) - 1).*tbc(nodes(b,3) > 0,6);
    surroundingconsidered(nodes(b,3) > 0,7) = (offset(nodes(b,3) > 0,1) - (x_dim)*(y_dim) - x_dim).*tbc(nodes(b,3) > 0,7);
    surroundingconsidered((nodes(b,3) > 0) & (nodes(b,2) > 0),8) = (offset((nodes(b,3) > 0) & (nodes(b,2) > 0),1) - (x_dim)*(y_dim) - (x_dim) -1).*tbc((nodes(b,3) > 0) & (nodes(b,2) > 0),8);
    %clear surrounding;
    clear tbc;

    tbcsum(tbcsum == 8) = 0;
    tbcsum(tbcsum == 4) = 0;
    tbcsum((distsum>2) & (mins > 1)) = 0;

        clear distsum;
    % get the surrounding elements which are to be considered for the shift


    % use the dummy centroid to make computations easier
    surroundingconsidered(surroundingconsidered == 0) = size(centroids,1);

    % calculate the combination of the centroids which are to be considered
    % for the shift
    centroidcomb = zeros(size(nodes));
    centroidcomb(:,1) = sum(reshape(centroids(surroundingconsidered,1), [], 8),2);
    centroidcomb(:,2) = sum(reshape(centroids(surroundingconsidered,2), [], 8),2);
    centroidcomb(:,3) = sum(reshape(centroids(surroundingconsidered,3), [], 8),2);
    clear surroundingconsidered;
    centroidcomb(tbcsum ~= 0,1) = centroidcomb(tbcsum ~=0,1)./tbcsum(tbcsum ~=0);
    centroidcomb(tbcsum ~= 0,2) = centroidcomb(tbcsum ~=0,2)./tbcsum(tbcsum ~=0);
    centroidcomb(tbcsum ~= 0,3) = centroidcomb(tbcsum ~=0,3)./tbcsum(tbcsum ~=0);

    % finally apply the shift
    nodes(tbcsum == 0,:) = points(tbcsum == 0,:);
    nodes(tbcsum ~= 0,:) = (1-sh)*nodes(tbcsum ~= 0,:) + sh*centroidcomb(tbcsum ~= 0,:);

function segmentation = ea_ft_datatype_segmentation(segmentation, varargin)

% FT_DATATYPE_SEGMENTATION describes the FieldTrip MATLAB structure for segmented
% voxel-based data and atlases. A segmentation can either be indexed or probabilistic
% (see below).
%
% A segmentation is a volumetric description which is usually derived from an anatomical
% MRI, which describes for each voxel the tissue type. It for example distinguishes
% between white matter, grey matter, csf, skull and skin. It is mainly used for masking
% in visualization, construction of volume conduction models and for construction of
% cortical sheets. An volume-based atlas is basically a very detailed segmentation with
% an anatomical label for each voxel.
%
% For example, the AFNI TTatlas+tlrc segmented brain atlas (which can be created
% with FT_READ_ATLAS) looks like this
%
%              dim: [161 191 141]        the size of the 3D volume in voxels
%        transform: [4x4 double]         affine transformation matrix for mapping the voxel coordinates to head coordinate system
%         coordsys: 'tal'                the transformation matrix maps the voxels into this (head) coordinate system
%             unit: 'mm'                 the units in which the coordinate system is expressed
%           brick0: [161x191x141 uint8]  integer values from 1 to N, the value 0 means unknown
%           brick1: [161x191x141 uint8]  integer values from 1 to M, the value 0 means unknown
%      brick0label: {Nx1 cell}
%      brick1label: {Mx1 cell}
%
% An example of a whole-brain anatomical MRI that was segmented using FT_VOLUMESEGMENT
% looks like this
%
%         dim: [256 256 256]         the size of the 3D volume in voxels
%   transform: [4x4 double]          affine transformation matrix for mapping the voxel coordinates to head coordinate system
%    coordsys: 'ctf'                 the transformation matrix maps the voxels into this (head) coordinate system
%        unit: 'mm'                  the units in which the coordinate system is expressed
%        gray: [256x256x256 double]  probabilistic map of the gray matter
%       white: [256x256x256 double]  probabilistic map of the white matter
%         csf: [256x256x256 double]  probabilistic map of the cerebrospinal fluid
%
% An example segmentation with binary values that can be used for construction of a
% BEM volume conduction model of the head looks like this
%
%           dim: [256 256 256]         the dimensionality of the 3D volume
%     transform: [4x4 double]          affine transformation matrix for mapping the voxel coordinates to head coordinate system
%      coordsys: 'ctf'                 the transformation matrix maps the voxels into this (head) coordinate system
%          unit: 'mm'                  the units in which the coordinate system is expressed
%         brain: [256x256x256 logical] binary map representing the voxels which belong to the brain
%         scalp: [256x256x256 logical] binary map representing the voxels which belong to the scalp
%         skull: [256x256x256 logical] binary map representing the voxels which belong to the skull
%
% The examples above demonstrate that a segmentation can be indexed, i.e. consisting of
% subsequent integer numbers (1, 2, ...) or probabilistic, consisting of real numbers
% ranging from 0 to 1 that represent probabilities between 0% and 100%. An extreme case
% is one where the probability is either 0 or 1, in which case the probability can be
% represented as a binary or logical array.
%
% The only difference to the volume data representation is that the segmentation
% structure contains the additional fields xxx and xxxlabel. See FT_DATATYPE_VOLUME for
% further details.
%
% Required fields:
%   - dim, transform
%
% Optional fields:
%   - coordsys, unit
%
% Deprecated fields:
%   - none
%
% Obsoleted fields:
%   - none
%
% Revision history:
% (2012/latest) The explicit distunction between the indexed and probabilistic
% representation was made. For the indexed representation the additional
% xxxlabel cell-array was introduced.
%
% (2005) The initial version was defined.
%
% See also FT_DATATYPE, FT_DATATYPE_VOLUME, FT_DATATYPE_PARCELLATION

% Copyright (C) 2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_datatype_segmentation.m 10212 2015-02-11 16:47:05Z roboos $

% get the optional input arguments, which should be specified as key-value pairs
version           = ft_getopt(varargin, 'version', 'latest');
segmentationstyle = ft_getopt(varargin, 'segmentationstyle');  % can be indexed or probabilistic
hasbrain          = ft_getopt(varargin, 'hasbrain', 'no');     % no means that it is not required, if present it won't be removed

% convert from string into boolean
hasbrain = istrue(hasbrain);

if strcmp(version, 'latest')
  segversion = '2012';
  volversion = 'latest';
  clear version
else
  segversion = version;
  volversion = version;
  clear version
end

if isempty(segmentation)
  return;
end

switch segversion
  case '2012'
    % determine whether the style of the input fields is probabilistic or indexed
    fn = fieldnames(segmentation);
    fn = setdiff(fn, 'inside'); % exclude the inside field from any conversions
    [indexed, probabilistic] = ea_determine_segmentationstyle(segmentation, fn, segmentation.dim);

    % ignore the fields that do not contain a segmentation
    sel = indexed | probabilistic;
    fn            = fn(sel);
    indexed       = indexed(sel);
    probabilistic = probabilistic(sel);

    % convert from an exclusive to cumulative representation
    % this is only only for demonstration purposes
    % for i=1:length(sel)
    %   segmentation.(fn{sel(i)}) = volumefillholes(segmentation.(fn{sel(i)}));
    % end

    [dum, i] = intersect(fn, {'scalp', 'skull', 'brain'});
    if numel(i)==3
      % put them in the preferred order
      fn(i) = {'scalp', 'skull', 'brain'};
    end
    [dum, i] = intersect(fn, {'skin', 'skull', 'brain'}); % this is not likely
    if numel(i)==3
      % put them in the preferred order
      fn(i) = {'skin', 'skull', 'brain'};
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ensure that the segmentation is internally consistent
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if any(probabilistic)
      segmentation = ea_fixsegmentation(segmentation, fn(probabilistic), 'probabilistic');
    end
    if any(indexed)
      segmentation = ea_fixsegmentation(segmentation, fn(indexed), 'indexed');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert the segmentation to the desired style
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(segmentationstyle)
      % keep it as it is
    elseif strcmp(segmentationstyle, 'indexed') && any(probabilistic)
      segmentation  = convert_segmentationstyle(segmentation, fn(probabilistic), segmentation.dim, 'indexed');
      indexed(probabilistic)       = true;  % these are now indexed
      probabilistic(probabilistic) = false; % these are now indexed
    elseif strcmp(segmentationstyle, 'probabilistic') && any(indexed)
      segmentation  = convert_segmentationstyle(segmentation, fn(indexed), segmentation.dim, 'probabilistic');
      probabilistic(indexed) = true;  % these are now probabilistic
      indexed(indexed)       = false; % these are now probabilistic
    end % converting between probabilistic and indexed

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add the brain if requested
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if hasbrain
      if all(indexed)
        fn = fieldnames(segmentation);
        sel = false(size(fn));
        for i=1:numel(fn)
          sel(i) = any(strcmp(fn, [fn{i} 'label']));
        end
        fn = fn(sel);

        if numel(fn)>1
          error('cannot construct a brain mask on the fly; this requires a single indexed representation');
        else
          seg      = segmentation.(fn{1});
          seglabel = segmentation.([fn{1} 'label']);
          if ~any(strcmp(seglabel, 'brain'))
            threshold = 0.5;
            smooth    = 5;
            % ensure that the segmentation contains the brain mask, if not then construct it from gray+white+csf
            if length(intersect(seglabel, {'gray' 'white' 'csf'}))~=3
              error('cannot construct a brain mask on the fly; this requires gray, white and csf');
            end
            gray  = seg==find(strcmp(seglabel, 'gray'));
            white = seg==find(strcmp(seglabel, 'white'));
            csf   = seg==find(strcmp(seglabel, 'csf'));
            brain = gray + white + csf;
            clear gray white csf seg
            brain = volumesmooth(brain,    smooth,    'brain');
            brain = volumethreshold(brain, threshold, 'brain');
            % store it in the output
            segmentation.brain = brain;
          end % try to construct the brain
        end

      elseif all(probabilistic)
        if ~isfield(segmentation, 'brain')
          if ~all(isfield(segmentation, {'gray' 'white' 'csf'}))
            error('cannot construct a brain mask on the fly; this requires gray, white and csf');
          end
          threshold = 0.5;
          smooth    = 5;
          % ensure that the segmentation contains the brain mask, if not then construct it from gray+white+csf tissue probability maps
          gray  = segmentation.gray;
          white = segmentation.white;
          csf   = segmentation.csf;
          brain = gray + white + csf;
          clear gray white csf
          brain = volumesmooth(brain,    smooth,    'brain');
          brain = volumethreshold(brain, threshold, 'brain');
          % store it in the output
          segmentation.brain = brain;
        end
      else
        error('cannot construct a brain mask on the fly; this requires a uniquely indexed or a uniquely probabilitic representation');
      end
    end % if hasbrain

  case '2005'
    % the only difference is that the indexed representation for xxx did not have the xxxlabel field prior to the 2012 version
    fn = fieldnames(segmentation);
    sel = ~cellfun(@isempty, regexp(fn, 'label$'));
    segmentation = rmfield(segmentation, fn(sel));
    % furthermore it corresponds to the oldest version of the volume representation
    volversion = '2003';

  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for segmentation datatype', segversion);
end

% the segmentation is a speciat type of volume structure, so ensure that it also fulfills the requirements for that
segmentation = ft_datatype_volume(segmentation, 'version', volversion);

function segmentation = ea_fixsegmentation(segmentation, fn, style)

% FIXSEGMENTATION is a helper function that ensures the segmentation to be internally
% consistent. It is used by ft_datatype_segmentation and ft_datatype_parcellation.
%
% % See also CONVERT_SEGMENTATIONSTYLE, DETERMINE_SEGMENTATIONSTYLE

switch style
  case 'indexed'

    for i=1:length(fn)
      indexval = unique(segmentation.(fn{i})(:));  % find the unique tissue types
      indexval = indexval(indexval~=0);            % these are the only ones that matter

      if any(indexval<0)
        error('an indexed representation cannot contain negative numbers');
      end

      if ~isfield(segmentation, [fn{i} 'label'])
        % ensure that the tissues have labels
        indexlabel = {};
        for j=1:length(indexval)
          indexlabel{indexval(j)} = sprintf('tissue %d', indexval(j));
        end
        segmentation.([fn{i} 'label']) = indexlabel;
      else
        % ensure that the tissue labels are consistent with the index values
        indexlabel = segmentation.([fn{i} 'label']);
        if numel(indexval)>numel(indexlabel)
          error('each index value should have a corresponding entry in %s', [fn{i} 'label']);
        elseif any(cellfun(@isempty, indexlabel(indexval)))
          error('each index value should have a corresponding entry in %s', [fn{i} 'label']);
        end
        % the checks above allow for the situation where
        %   indexval   = [1 2 4]
        %   indexlabel = {'a', 'b', 'c', 'd'} or {'a', 'b', [], 'd'}
        % which happens if the segmentation unexpectedly does not contain a certain tissue type
      end

      % ensure that the indices are subsequent integers, i.e. [1 2 3] rather than [1 2 4]
      for j=1:length(indexval)
        tmp = segmentation.(fn{i});
        tmp(tmp==indexval(j)) = j;
        segmentation.(fn{i}) = tmp;
      end
      segmentation.([fn{i} 'label']) = segmentation.([fn{i} 'label'])(indexval);
    end
    clear tmp indexval indexlabel

  case 'probabilistic'

    % convert from a cumulative to an exclusive representation
    within = false(length(fn));
    if length(fn)>4
      % test for each tissue whether it is overlapping with or contained in each other tissue
      warning('more than 4 tissue types, this may take a while');
    end
    for i=1:length(fn)
      segi = segmentation.(fn{i})>0;
      for j=1:length(fn)
        if i==j
          % don't test for self-overlap
          continue
        end
        if ~any(segi(:))
          % don't bother to test completely empty segmentations
          continue
        end
        segj = segmentation.(fn{j})>0;
        within(i,j) = all(segj(segi(:))); % segi is fully contained in segj
        if i~=j && within(i,j)
          fprintf('the %s is fully contained in the %s, removing it from the %s\n', fn{i}, fn{j}, fn{j});
          segmentation.(fn{j})(segi) = 0;
        end
      end
    end
    clear segi segj within

  otherwise
    error('unsupported style "%s"', style);
end

function ea_ft_postamble(cmd, varargin)

% FT_POSTAMBLE is a helper function that is included in many of the FieldTrip
% functions and which takes care of some general settings and operations at the end
% of the function.
%
% This ft_postamble m-file is a function, but internally it executes a number of
% private scripts in the callers workspace. This allows the private script to access
% the variables in the callers workspace and behave as if the script were included as
% a header file in C-code.
%
% See also FT_PREAMBLE

% Copyright (C) 2011-2012, Robert Oostenveld, DCCN
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_postamble.m 9573 2014-05-21 15:17:03Z roboos $

% ideally this would be a script, because the local variables would then be
% shared with the calling function. Instead, this is a function which then
% passes the variables explicitely to another script which is eval'ed.

% the following section ensures that these scripts are included as
% dependencies when using the MATLAB compiler
%
%#function ft_postamble_debug
%#function ft_postamble_trackconfig
%#function ft_postamble_provenance
%#function ft_postamble_previous
%#function ft_postamble_history
%#function ft_postamble_savevar
%#function ft_postamble_randomseed

global ft_default

% this is a trick to pass the input arguments into the ft_postamble_xxx script
ft_default.postamble = varargin;

if exist(['ft_postamble_' cmd], 'file')
  evalin('caller', ['ft_postamble_' cmd]);
end

if isfield(ft_default, 'postamble')
  % the postamble field should not remain in the ft_default structure
  ft_default = rmfield(ft_default, 'postamble');
end

function val = ea_ft_getopt(opt, key, default, emptymeaningful)

% FT_GETOPT gets the value of a specified option from a configuration structure
% or from a cell-array with key-value pairs.
%
% Use as
%   val = ft_getopt(s, key, default, emptymeaningful)
% where the input values are
%   s               = structure or cell-array
%   key             = string
%   default         = any valid MATLAB data type
%   emptymeaningful = boolean value (optional, default = 0)
%
% If the key is present as field in the structure, or as key-value
% pair in the cell-array, the corresponding value will be returned.
%
% If the key is not present, ft_getopt will return an empty array.
%
% If the key is present but has an empty value, then the emptymeaningful
% flag specifies whether the empty value or the default value should
% be returned. If emptymeaningful==true, then an empty array will be
% returned. If emptymeaningful==false, then the specified default will
% be returned.
%
% See also FT_SETOPT, FT_CHECKOPT

% Copyright (C) 2011-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_getopt.m 7123 2012-12-06 21:21:38Z roboos $

if nargin<3
  default = [];
end

if nargin < 4
  emptymeaningful = 0;
end

if isa(opt, 'struct') || isa(opt, 'config')
  % get the key-value from the structure
  fn = fieldnames(opt);
  if ~any(strcmp(key, fn))
    val = default;
  else
    val = opt.(key);
  end

elseif isa(opt, 'cell')
  % get the key-value from the cell-array
  if mod(length(opt),2)
    error('optional input arguments should come in key-value pairs, i.e. there should be an even number');
  end

  % the 1st, 3rd, etc. contain the keys, the 2nd, 4th, etc. contain the values
  keys = opt(1:2:end);
  vals = opt(2:2:end);

  % the following may be faster than cellfun(@ischar, keys)
  valid = false(size(keys));
  for i=1:numel(keys)
    valid(i) = ischar(keys{i});
  end

  if ~all(valid)
    error('optional input arguments should come in key-value pairs, the optional input argument %d is invalid (should be a string)', i);
  end

  hit = find(strcmpi(key, keys));
  if isempty(hit)
    % the requested key was not found
    val = default;
  elseif length(hit)==1
    % the requested key was found
    val = vals{hit};
  else
    error('multiple input arguments with the same name');
  end

elseif isempty(opt)
  % no options are specified, return default
  val = default;
end % isstruct or iscell or isempty

if isempty(val) && ~isempty(default) && ~emptymeaningful
  % use the default value instead of the empty input that was specified:
  % this applies for example if you do functionname('key', []), where
  % the empty is meant to indicate that the user does not know or care
  % what the value is
  val = default;
end

function [pnt, tri] = ea_makesphere(numvertices)

if isempty(numvertices)
  [pnt,tri] = icosahedron162;
  fprintf('using the mesh specified by icosaedron162\n');
elseif numvertices==42
  [pnt,tri] = icosahedron42;
  fprintf('using the mesh specified by icosaedron%d\n',size(pnt,1));
elseif numvertices==162
  [pnt,tri] = icosahedron162;
  fprintf('using the mesh specified by icosaedron%d\n',size(pnt,1));
elseif numvertices==642
  [pnt,tri] = icosahedron642;
  fprintf('using the mesh specified by icosaedron%d\n',size(pnt,1));
elseif numvertices==2562
  [pnt,tri] = icosahedron2562;
  fprintf('using the mesh specified by icosaedron%d\n',size(pnt,1));
else
  [pnt, tri] = msphere(numvertices);
  fprintf('using the mesh specified by msphere with %d vertices\n',size(pnt,1));
end

function ea_ft_preamble(cmd, varargin)

% FT_PREAMBLE is a helper function that is included in many of the FieldTrip
% functions and which takes care of some general settings and operations at the
% begin of the function.
%
% This ft_preamble m-file is a function, but internally it executes a
% number of private scripts in the callers workspace. This allows the
% private script to access the variables in the callers workspace and
% behave as if the script were included as a header file in C-code.
%
% See also FT_POSTAMBLE

% Copyright (C) 2011-2012, Robert Oostenveld, DCCN
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_preamble.m 9520 2014-05-14 09:33:28Z roboos $

% ideally this would be a script, because the local variables would then be
% shared with the calling function. Instead, this is a function which then
% passes the variables explicitely to another script which is eval'ed.

% the following section ensures that these scripts are included as
% dependencies when using the MATLAB compiler
%
%#function ft_preamble_init
%#function ft_preamble_debug
%#function ft_preamble_trackconfig
%#function ft_preamble_provenance
%#function ft_preamble_loadvar
%#function ft_preamble_randomseed

global ft_default

% this is a trick to pass the input arguments into the ft_preamble_xxx script
ft_default.preamble = varargin;

if exist(['ft_preamble_' cmd], 'file')
  evalin('caller', ['ft_preamble_' cmd]);
end

if isfield(ft_default, 'preamble')
  % the preamble field should not remain in the ft_default structure
  ft_default = rmfield(ft_default, 'preamble');
end

function ea_ft_defaults

% FT_DEFAULTS (ending with "s") sets some general settings in the global variable
% ft_default (without the "s") and takes care of the required path settings. This
% function is called at the begin of all FieldTrip functions.
%
% The configuration defaults are stored in the global "ft_default" structure.
% The ft_checkconfig function that is called by many FieldTrip functions will
% merge this global ft_default structure with the cfg ctructure that you pass to
% the FieldTrip function that you are calling.
%
% The global options and their default values are
%   ft_default.trackconfig    string, can be cleanup, report, off (default = 'off')
%   ft_default.checkconfig    string, can be pedantic, loose, silent (default = 'loose')
%   ft_default.checksize      number in bytes, can be inf (default = 1e5)
%   ft_default.showcallinfo   string, can be yes or no (default = 'yes')
%   ft_default.debug          string, can be 'display', 'displayonerror', 'displayonsuccess',
%                             'save', 'saveonerror', saveonsuccess' or 'no' (default = 'no')
%
% See also FT_HASTOOLBOX, FT_CHECKCONFIG

% Note that this should be a function and not a script, otherwise the
% ft_hastoolbox function appears not be found in fieldtrip/private.

% Copyright (C) 2009-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_defaults.m 9674 2014-06-24 07:49:49Z eelspa $

global ft_default
persistent initialized

% Set the defaults in a global variable, ft_checkconfig will copy these over into the local configuration.
% Note that ft_getopt might not be available on the path at this moment and can therefore not yet be used.

if ~isfield(ft_default, 'trackconfig'),    ft_default.trackconfig    = 'off';    end % cleanup, report, off
if ~isfield(ft_default, 'checkconfig'),    ft_default.checkconfig    = 'loose';  end % pedantic, loose, silent
if ~isfield(ft_default, 'checksize'),      ft_default.checksize      = 1e5;      end % number in bytes, can be inf
if ~isfield(ft_default, 'showcallinfo'),   ft_default.showcallinfo   = 'yes';    end % yes or no, this is used in ft_pre/postamble_provenance
if ~isfield(ft_default, 'debug'),          ft_default.debug          = 'no';     end % no, save, saveonerror, display, displayonerror, this is used in ft_pre/postamble_debug

% these options allow to disable parts of the provenance
if ~isfield(ft_default, 'trackcallinfo'),  ft_default.trackcallinfo  = 'yes';    end % yes or no
if ~isfield(ft_default, 'trackdatainfo'),  ft_default.trackdatainfo  = 'no';     end % yes or no, this is still under development
if ~isfield(ft_default, 'trackparaminfo'), ft_default.trackparaminfo = 'no';     end % yes or no, this is still under development

% track whether we have executed ft_defaults already. Note that we should
% not use ft_default itself directly, because the user might have set stuff
% in that struct already before ft_defaults is called for the first time.
if ~isempty(initialized) && exist('ft_hastoolbox', 'file')
  return;
end

% Ensure that the path containing ft_defaults is on the path.
% This allows people to do "cd path_to_fieldtrip; ft_defaults"
ftPath = fileparts(mfilename('fullpath')); % get path, strip away 'ft_defaults'
ftPath = strrep(ftPath, '\', '\\');
if isempty(regexp(path, [ftPath pathsep '|' ftPath '$'], 'once'))
  warning('FieldTrip is not yet on your MATLAB path, adding %s', strrep(ftPath, '\\', '\'));
  addpath(ftPath);
end

if ~isdeployed

  % Some people mess up their path settings and then have
  % different versions of certain toolboxes on the path.
  % The following will issue a warning
  ea_checkMultipleToolbox('FieldTrip',           'ft_defaults.m');
  ea_checkMultipleToolbox('spm',                 'spm.m');
  ea_checkMultipleToolbox('mne',                 'fiff_copy_tree.m');
  ea_checkMultipleToolbox('eeglab',              'eeglab2fieldtrip.m');
  ea_checkMultipleToolbox('dipoli',              'write_tri.m');
  ea_checkMultipleToolbox('eeprobe',             'read_eep_avr.mexa64');
  ea_checkMultipleToolbox('yokogawa',            'GetMeg160ChannelInfoM.p');
  ea_checkMultipleToolbox('simbio',              'sb_compile_vista.m');
  ea_checkMultipleToolbox('fns',                 'fns_region_read.m');
  ea_checkMultipleToolbox('bemcp',               'bem_Cii_cst.mexa64');
  ea_checkMultipleToolbox('bci2000',             'load_bcidat.m');
  ea_checkMultipleToolbox('openmeeg',            'openmeeg_helper.m');
  ea_checkMultipleToolbox('freesurfer',          'vox2ras_ksolve.m');
  ea_checkMultipleToolbox('fastica',             'fastica.m');
  ea_checkMultipleToolbox('besa',                'readBESAmul.m');
  ea_checkMultipleToolbox('neuroshare',          'ns_GetAnalogData.m');
  ea_checkMultipleToolbox('ctf',                 'setCTFDataBalance.m');
  ea_checkMultipleToolbox('afni',                'WriteBrikHEAD.m');
  ea_checkMultipleToolbox('gifti',               '@gifti/display.m');
  ea_checkMultipleToolbox('sqdproject',          'sqdread.m');
  ea_checkMultipleToolbox('xml4mat',             'xml2mat.m');
  ea_checkMultipleToolbox('cca',                 'ccabss.m');
  ea_checkMultipleToolbox('bsmart',              'armorf.m');
  ea_checkMultipleToolbox('iso2mesh',            'iso2meshver.m');
  ea_checkMultipleToolbox('bct',                 'degrees_und.m');
  ea_checkMultipleToolbox('yokogawa_meg_reader', 'getYkgwHdrEvent.p');
  ea_checkMultipleToolbox('biosig',              'sopen.m');
  ea_checkMultipleToolbox('icasso',              'icassoEst.m');

  if isempty(which('ft_hastoolbox'))
    % the fieldtrip/utilities directory contains the ft_hastoolbox function
    % which is required for the remainder of this script
    addpath(fullfile(fileparts(which('ft_defaults')), 'utilities'));
  end

  try
    % external/signal directory contains alternative implementations of some signal processing functions
    addpath(fullfile(fileparts(which('ft_defaults')), 'external', 'signal'));
  end

  try
    % some alternative implementations of statistics functions
    addpath(fullfile(fileparts(which('ft_defaults')), 'external', 'stats'));
  end

  try
    % this directory contains various functions that were obtained from elsewere, e.g. MATLAB file exchange
    ea_ft_hastoolbox('fileexchange', 3, 1); % not required
  end

  try
    % this directory contains the backward compatibility wrappers for the ft_xxx function name change
    ea_ft_hastoolbox('compat', 3, 1); % not required
  end

  try
    % these directories contain functions that were added to MATLAB in
    % recent versions to replace an older function.
    if matlabversion(-Inf, '2011b')
      ea_ft_hastoolbox('compat/matlablt2012a', 2, 1);
    end
  end

  try
    % these contains template layouts, neighbour structures, MRIs and cortical meshes
    ea_ft_hastoolbox('template/layout', 1, 1);
    ea_ft_hastoolbox('template/anatomy', 1, 1);
    ea_ft_hastoolbox('template/headmodel', 1, 1);
    ea_ft_hastoolbox('template/electrode', 1, 1);
    ea_ft_hastoolbox('template/neighbours', 1, 1);
    ea_ft_hastoolbox('template/sourcemodel', 1, 1);
  end

  try
    % this is used in statistics
    ea_ft_hastoolbox('statfun', 1, 1);
  end

  try
    % this is used in definetrial
    ea_ft_hastoolbox('trialfun', 1, 1);
  end

  try
    % this contains the low-level reading functions
    ea_ft_hastoolbox('fileio', 1, 1);
  end

  try
    % this is for filtering time-series data
    ea_ft_hastoolbox('preproc', 1, 1);
  end

  try
    % this contains forward models for the EEG and MEG volume conduction problem
    ea_ft_hastoolbox('forward', 1, 1);
  end

  try
    % numerous functions depend on this module
    ea_ft_hastoolbox('inverse', 1, 1);
  end

  try
    % this contains intermediate-level plotting functions, e.g. multiplots and 3-d objects
    ea_ft_hastoolbox('plotting', 1, 1);
  end

  try
    % this contains the functions to compute connecitivy metrics
    ea_ft_hastoolbox('connectivity', 1,1);
  end

  try
    % this contains the functions for spike and spike-field analysis
    ea_ft_hastoolbox('spike', 1,1);
  end

  try
    % this contains specific code and examples for realtime processing
    ea_ft_hastoolbox('realtime/example', 3, 1);    % not required
    ea_ft_hastoolbox('realtime/online_mri', 3, 1); % not required
    ea_ft_hastoolbox('realtime/online_meg', 3, 1); % not required
    ea_ft_hastoolbox('realtime/online_eeg', 3, 1); % not required
  end

  try
    % this contains intermediate-level functions for spectral analysis
    ea_ft_hastoolbox('specest', 1, 1);
  end

end

% remember that the function has executed in a persistent variable
initialized = true;

function ea_checkMultipleToolbox(toolbox, keyfile)
% do nothing (delete this function one day)..

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

function [lf] = ea_leadfield_simbio(pos, vol)

% leadfield_simbio leadfields for a set of dipoles
%
% [lf] = leadfield_simbio(pos, vol);
%
% with input arguments
%   pos     a matrix of dipole positions
%           there can be 'deep electrodes' too!
%   vol     contains a FE volume conductor (output of ft_prepare_vol_sens)
%
% the output lf is the leadfield matrix of dimensions m (rows) x n*3 (cols)

% copyright (c) 2012, Johannes Vorwerk

try
    lf = zeros(size(3*pos,1),size(vol.transfer,1));
    dir = diag([1,1,1]);
    for i=1:size(pos,1)
        locpos = repmat(pos(i,:),3,1);
        rhs = ea_sb_rhs_venant(locpos,dir,vol);
        lf((3*(i-1)+1):(3*(i-1)+3),:) = (vol.transfer * rhs)';
    end
    lf = lf';
catch
  warning('an error occurred while running simbio');
  rethrow(lasterror)
end

function rhs = ea_sb_rhs_venant(pos,dir,vol)

% SB_RHS_VENANT
%
% $Id: sb_rhs_venant.m 8776 2013-11-14 09:04:48Z roboos $

%find node closest to source position
next_nd = ea_sb_get_next_nd(pos,vol.pos);
%find nodes neighbouring closest node
if isfield(vol,'tet')
    ven_nd = ea_sb_get_ven_nd(next_nd,vol.tet);
elseif isfield(vol,'hex')
    ven_nd = ea_sb_get_ven_nd(next_nd,vol.hex);
else
    error('No connectivity information given!');
end
%calculate rhs matrix
loads = ea_sb_calc_ven_loads(pos,dir,ven_nd,vol.pos);
%assign values in sparse matrix
i = reshape(ven_nd',[],1);
j = reshape(repmat([1:size(ven_nd,1)],size(ven_nd,2),1),[],1);
loads = reshape(loads',[],1);
j = j(i~=0);
loads = loads(i~=0);
i = i(i~=0);
rhs = sparse(i,j,loads,size(vol.pos,1),size(pos,1));


function next_nd = ea_sb_get_next_nd(pos,node);
next_nd = zeros(size(pos,1),1);
for i=1:size(pos,1)
    [dist, next_nd(i)] = min(sum(bsxfun(@minus,node,pos(i,:)).^2,2));
end

function loads = ea_sb_calc_ven_loads(pos,dir,ven_nd,node);

% SB_CALC_VEN_LOADS
%
% $Id: sb_calc_ven_loads.m 8776 2013-11-14 09:04:48Z roboos $

%aref setzen
aref = 20;
%lambda setzen
lambda = 10e-6;
%r setzen
r = 1;

loads = zeros(size(pos,1),size(ven_nd,2));
for i=1:size(pos,1);
    x = bsxfun(@minus,node(ven_nd(i,ven_nd(i,:)~=0),:),pos(i,:))./aref;
    X = zeros(9,size(x,1));
    X(1:3:7,:) = ones(3,size(x,1));
    X(2:3:8,:) = x';
    X(3:3:9,:) = (x.^2)';
    T = zeros(9,1);
    T(1:3:7) = 0;
    T(2:3:8) = dir(i,:)./aref;
    T(3:3:9) = 0;
    W = zeros(size(x,1)*3,size(x,1));
    W(1:size(x,1),:) = diag(x(:,1).^r);
    W(size(x,1)+1:2*size(x,1),:) = diag(x(:,2).^r);
    W(size(x,1)*2+1:3*size(x,1),:) = diag(x(:,3).^r);
    tmp = ((X')*X + lambda*(W')*W) \ X'*T;
    loads(i,ven_nd(i,:)~=0) = tmp;
end


function ven_nd = ea_sb_get_ven_nd(next_nd,elem);
ven_nd = zeros(size(next_nd,1),1);
for i=1:size(next_nd,1)
    [tmp1,tmp2] = find(elem == next_nd(i));
    tmp = unique(elem(tmp1,:));
    %tmp = tmp(tmp~=next_nd(i)); seems like this is not done in the
    %original.
    if(length(tmp) > size(ven_nd,2))
        ven_nd = [ven_nd, zeros(size(ven_nd,1),length(tmp)-size(ven_nd,2))];
    end
    ven_nd(i,1:length(tmp)) = tmp;
end



function [output] = ea_ft_transform_geometry(transform, input)

% FT_TRANSFORM_GEOMETRY applies a homogeneous coordinate transformation to
% a structure with geometric information. These objects include:
%  - volume conductor geometry, consisting of a mesh, a set of meshes, a
%      single sphere, or multiple spheres.
%  - gradiometer of electrode structure containing sensor positions and
%      coil orientations (for MEG).
%  - headshape description containing positions in 3D space.
%  - sourcemodel description containing positions and optional orientations
%      in 3D space.
%
% The units in which the transformation matrix is expressed are assumed to
% be the same units as the units in which the geometric object is
% expressed. Depending on the input object, the homogeneous transformation
% matrix should be limited to a rigid-body translation plus rotation
% (MEG-gradiometer array), or to a rigid-body translation plus rotation
% plus a global rescaling (volume conductor geometry).
%
% Use as
%   output = ft_transform_geometry(transform, input)

% Copyright (C) 2011, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_transform_geometry.m$

% flg rescaling check
allowscaling = ~ft_senstype(input, 'meg');

% determine the rotation matrix
rotation = eye(4);
rotation(1:3,1:3) = transform(1:3,1:3);

if any(abs(transform(4,:)-[0 0 0 1])>100*eps)
  error('invalid transformation matrix');
end

if ~allowscaling
  % allow for some numerical imprecision
  if abs(det(rotation)-1)>1e-6%100*eps
  %if abs(det(rotation)-1)>100*eps  % allow for some numerical imprecision
    error('only a rigid body transformation without rescaling is allowed');
  end
end

if allowscaling
  % FIXME build in a check for uniform rescaling probably do svd or so
  % FIXME insert check for nonuniform scaling, should give an error
end

tfields   = {'pos' 'pnt' 'o' 'chanpos' 'chanposorg' 'coilpos' 'elecpos', 'nas', 'lpa', 'rpa', 'zpoint'}; % apply rotation plus translation
rfields   = {'ori' 'nrm' 'coilori'}; % only apply rotation
mfields   = {'transform'};           % plain matrix multiplication
recfields = {'fid' 'bnd' 'orig'};    % recurse into these fields
% the field 'r' is not included here, because it applies to a volume
% conductor model, and scaling is not allowed, so r will not change.

fnames    = fieldnames(input);
for k = 1:numel(fnames)
    if ~isempty(input.(fnames{k}))
        if any(strcmp(fnames{k}, tfields))
            input.(fnames{k}) = ea_apply(transform, input.(fnames{k}));
        elseif any(strcmp(fnames{k}, rfields))
            input.(fnames{k}) = ea_apply(rotation, input.(fnames{k}));
        elseif any(strcmp(fnames{k}, mfields))
            input.(fnames{k}) = transform*input.(fnames{k});
        elseif any(strcmp(fnames{k}, recfields))
            for j = 1:numel(input.(fnames{k}))
                input.(fnames{k})(j) = ft_transform_geometry(transform, input.(fnames{k})(j));
            end
        else
            % do nothing
        end
    end
end
output = input;
return;

function [new] = ea_apply(transform, old)
old(:,4) = 1;
new = old * transform';
new = new(:,1:3);

function [type] = ea_ft_senstype(input, desired)

% FT_SENSTYPE determines the type of acquisition device by looking at the
% channel names and comparing them with predefined lists.
%
% Use as
%   [type] = ft_senstype(sens)
%   [flag] = ft_senstype(sens, desired)
%
% The output type can be any of the following
%   'ctf64'
%   'ctf151'
%   'ctf151_planar'
%   'ctf275'
%   'ctf275_planar'
%   'bti148'
%   'bti148_planar'
%   'bti248'
%   'bti248_planar'
%   'bti248grad'
%   'bti248grad_planar'
%   'itab28'
%   'itab153'
%   'itab153_planar'
%   'yokogawa9'
%   'yokogawa64'
%   'yokogawa64_planar'
%   'yokogawa160'
%   'yokogawa160_planar'
%   'yokogawa440'
%   'yokogawa440'_planar
%   'ext1020' (this includes eeg1020, eeg1010 and eeg1005)
%   'neuromag122'
%   'neuromag306'
%   'babysquid74'
%   'egi32'
%   'egi64'
%   'egi128'
%   'egi256'
%   'biosemi64'
%   'biosemi128'
%   'biosemi256'
%   'neuralynx'
%   'plexon'
%   'artinis'
%   'eeg' (this was called 'electrode' in older versions)
%   'meg' (this was called 'magnetometer' in older versions)
%   'nirs'
%
% The optional input argument for the desired type can be any of the above,
% or any of the following
%   'eeg'
%   'meg'
%   'meg_planar'
%   'meg_axial'
%   'ctf'
%   'bti'
%   'neuromag'
%   'yokogawa'
% If you specify the desired type, this function will return a boolean
% true/false depending on the input data.
%
% Besides specifiying a sensor definition (i.e. a grad or elec structure,
% see FT_DATATYPE_SENS), it is also possible to give a data structure
% containing a grad or elec field, or giving a list of channel names (as
% cell-arrray). So assuming that you have a FieldTrip data structure, any
% of the following calls would also be fine.
%   ft_senstype(hdr)
%   ft_senstype(data)
%   ft_senstype(data.label)
%   ft_senstype(data.grad)
%   ft_senstype(data.grad.label)
%
% See also FT_SENSLABEL, FT_CHANTYPE, FT_READ_SENS, FT_COMPUTE_LEADFIELD, FT_DATATYPE_SENS

% Copyright (C) 2007-2014, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_senstype.m 10340 2015-04-17 14:10:04Z jorhor $

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

% this is to avoid a recursion loop
persistent recursion
if isempty(recursion)
  recursion = false;
end

if iscell(input) && numel(input)<4 && ~all(cellfun(@ischar, input))
  % this might represent combined EEG, ECoG and/or MEG
  type = cell(size(input));
  if nargin<2
    desired = cell(size(input)); % empty elements
  end
  for i=1:numel(input)
    type{i} = ea_ft_senstype(input{i}, desired{i});
  end
  return
end

if nargin<2
  % ensure that all input arguments are defined
  desired = [];
end

current_argin = {input, desired};
if isequal(current_argin, previous_argin)
  % don't do the type detection again, but return the previous output from cache
  type = previous_argout{1};
  return
end

isdata   = isa(input, 'struct')  && (isfield(input, 'hdr') || isfield(input, 'time') || isfield(input, 'freq') || isfield(input, 'grad') || isfield(input, 'elec') || isfield(input, 'opto'));
isheader = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'Fs');
isgrad   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'pnt')  &&  isfield(input, 'ori'); % old style
iselec   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'pnt')  && ~isfield(input, 'ori'); % old style
isgrad   = (isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'coilpos')) || isgrad;             % new style
iselec   = (isa(input, 'struct') && isfield(input, 'label') && isfield(input, 'elecpos')) || iselec;             % new style
isnirs   = isa(input, 'struct')  && isfield(input, 'label') && isfield(input, 'transceiver');
islabel  = isa(input, 'cell')    && ~isempty(input) && isa(input{1}, 'char');
haslabel = isa(input, 'struct')  && isfield(input, 'label');

if ~(isdata || isheader || isgrad || iselec || isnirs || islabel || haslabel) && isfield(input, 'hdr')
  input    = input.hdr;
  isheader = true;
end

if isdata
  % the input may be a data structure which then contains a grad/elec structure, a header or only the labels
  % preferably look at the data and not the header for the grad, because it might be re-balanced and/or planar
  if isfield(input, 'grad')
    sens   = input.grad;
    isgrad = true;
  elseif isfield(input, 'elec')
    sens   = input.elec;
    iselec = true;
  elseif issubfield(input, 'hdr.grad')
    sens   = input.hdr.grad;
    isgrad = true;
  elseif issubfield(input, 'hdr.elec')
    sens   = input.hdr.elec;
    iselec = true;
  elseif issubfield(input, 'hdr.opto')
    sens   = input.hdr.opto;
    isnirs = true;
  elseif issubfield(input, 'hdr.label')
    sens.label = input.hdr.label;
    islabel    = true;
  elseif isfield(input, 'label')
    sens.label = input.label;
    islabel    = true;
  end

elseif isheader
  if isfield(input, 'grad')
    sens   = input.grad;
    isgrad = true;
  elseif isfield(input, 'elec')
    sens   = input.elec;
    iselec = true;
  elseif isfield(input, 'opto')
    sens   = input.opto;
    isnirs = true;
  elseif isfield(input, 'label')
    sens.label = input.label;
    islabel    = true;
  end

elseif isgrad
  sens = input;

elseif iselec
  sens = input;

elseif isnirs
  sens = input;

elseif islabel
  sens.label = input;

elseif haslabel
  % it does not resemble anything that we had expected at this location, but it does have channel labels
  % the channel labels can be used to determine the type of sensor array
  sens.label = input.label;
  islabel    = true;

else
  sens = [];
end


if isfield(sens, 'type')
  % preferably the structure specifies its own type
  type = sens.type;

  % do not make a distinction between the neuromag data with or without space in the channel names
  if strcmp(type, 'neuromag306alt')
    type = 'neuromag306';
  elseif strcmp(type, 'neuromag122alt')
    type = 'neuromag122';
  end

elseif isfield(input, 'nChans') && input.nChans==1 && isfield(input, 'label') && ~isempty(regexp(input.label{1}, '^csc', 'once'))
  % this is a single channel header that was read from a Neuralynx file, might be fcdc_matbin or neuralynx_nsc
  type = 'neuralynx';

elseif issubfield(input, 'orig.FileHeader') &&  issubfield(input, 'orig.VarHeader')
  % this is a complete header that was read from a Plexon *.nex file using read_plexon_nex
  type = 'plexon';

elseif issubfield(input, 'orig.stname')
  % this is a complete header that was read from an ITAB dataset
  type = 'itab';

elseif issubfield(input, 'orig.sys_name')
  % this is a complete header that was read from a Yokogawa dataset
  if strcmp(input.orig.sys_name, '9ch Biomagnetometer System') || input.orig.channel_count<20
    % this is the small animal system that is installed at the UCL Ear Institute
    % see http://www.ucl.ac.uk/news/news-articles/0907/09070101
    type = 'yokogawa9';
  elseif input.orig.channel_count<160
    type = 'yokogawa64';
  elseif input.orig.channel_count<300
    type = 'yokogawa160';
  else
    % FIXME this might fail if there are many bad channels
    type = 'yokogawa440';
  end

elseif issubfield(input, 'orig.FILE.Ext') && strcmp(input.orig.FILE.Ext, 'edf')
  % this is a complete header that was read from an EDF or EDF+ dataset
  type = 'eeg';

else
  % start with unknown, then try to determine the proper type by looking at the labels
  type = 'unknown';

  if isgrad && isfield(sens, 'type')
    type = sens.type;

  elseif isgrad
    % this looks like MEG
    % revert the component balancing that was previously applied
    if isfield(sens, 'balance') && strcmp(sens.balance.current, 'comp')
      sens = ea_undobalancing(sens);
    end

    % determine the type of magnetometer/gradiometer system based on the channel names alone
    % this uses a recursive call to the "islabel" section further down
    type = ea_ft_senstype(sens.label);
    %sens_temp.type = type;
    if strcmp(type, 'unknown') %|| ~ft_senstype(sens_temp,'meg')
      % although we don't know the type, we do know that it is MEG
      type = 'meg';
    end

  elseif iselec
    % this looks like EEG

    % determine the type of eeg/acquisition system based on the channel names alone
    % this uses a recursive call to the "islabel" section further down
    type = ea_ft_senstype(sens.label);
    %sens_temp.type = type;
    if strcmp(type, 'unknown') %|| ~ft_senstype(sens_temp,'eeg')
      % although we don't know the type, we do know that it is EEG
      type = 'eeg';
    end

  elseif isnirs
    % this looks like NIRS

    % determine the type of eeg/acquisition system based on the channel names alone
    % this uses a recursive call to the "islabel" section further down
    type = ft_senstype(sens.label);
    %sens_temp.type = type;
    if strcmp(type, 'unknown') %|| ~ft_senstype(sens_temp,'eeg')
      % although we don't know the type, we do know that it is EEG
      type = 'nirs';
    end

  elseif islabel
    % look only at the channel labels
    if     (mean(ismember(ea_ft_senslabel('ant128'),         sens.label)) > 0.8)
      type = 'ant128';
    elseif (mean(ismember(ea_ft_senslabel('ctf275'),         sens.label)) > 0.8)
      type = 'ctf275';
    elseif (mean(ismember(ea_ft_senslabel('ctfheadloc'),     sens.label)) > 0.8)  % look at the head localization channels
      type = 'ctf275';
    elseif (mean(ismember(ea_ft_senslabel('ctf151'),         sens.label)) > 0.8)
      type = 'ctf151';
    elseif (mean(ismember(ea_ft_senslabel('ctf64'),          sens.label)) > 0.8)
      type = 'ctf64';
    elseif (mean(ismember(ea_ft_senslabel('ctf275_planar'),  sens.label)) > 0.8)
      type = 'ctf275_planar';
    elseif (mean(ismember(ea_ft_senslabel('ctf151_planar'),  sens.label)) > 0.8)
      type = 'ctf151_planar';
    elseif (mean(ismember(ea_ft_senslabel('bti248'),         sens.label)) > 0.8) % note that it might also be a bti248grad system
      type = 'bti248';
    elseif (mean(ismember(ea_ft_senslabel('bti148'),         sens.label)) > 0.8)
      type = 'bti148';
    elseif (mean(ismember(ea_ft_senslabel('bti248_planar'),  sens.label)) > 0.8) % note that it might also be a bti248grad_planar system
      type = 'bti248_planar';
    elseif (mean(ismember(ea_ft_senslabel('bti148_planar'),  sens.label)) > 0.8)
      type = 'bti148_planar';
    elseif (mean(ismember(ea_ft_senslabel('itab28'),         sens.label)) > 0.8)
      type = 'itab28';
    elseif (mean(ismember(ea_ft_senslabel('itab153'),        sens.label)) > 0.8)
      type = 'itab153';
    elseif (mean(ismember(ea_ft_senslabel('itab153_planar'), sens.label)) > 0.8)
      type = 'itab153_planar';

      % the order is important for the different yokogawa systems, because they all share the same channel names
    elseif (mean(ismember(ea_ft_senslabel('yokogawa440'),        sens.label)) > 0.7)
      type = 'yokogawa440';
    elseif (mean(ismember(ea_ft_senslabel('yokogawa440_planar'), sens.label)) > 0.7)
      type = 'yokogawa440_planar';
    elseif (mean(ismember(ea_ft_senslabel('yokogawa160'),        sens.label)) > 0.4)
      type = 'yokogawa160';
    elseif (mean(ismember(ea_ft_senslabel('yokogawa160_planar'), sens.label)) > 0.4)
      type = 'yokogawa160_planar';
    elseif (mean(ismember(ea_ft_senslabel('yokogawa64'),         sens.label)) > 0.4)
      type = 'yokogawa64';
    elseif (mean(ismember(ea_ft_senslabel('yokogawa64_planar'),  sens.label)) > 0.4)
      type = 'yokogawa64_planar';
    elseif all(ismember(ea_ft_senslabel('yokogawa9'),            sens.label))
      type = 'yokogawa9';

    elseif any(mean(ismember(ea_ft_senslabel('neuromag306'),     sens.label)) > 0.4) % there are two possibilities for the channel labels: with and without a space
      type = 'neuromag306';
    elseif any(mean(ismember(ea_ft_senslabel('neuromag122'),     sens.label)) > 0.4) % there are two possibilities for the channel labels: with and without a space
      type = 'neuromag122';

    elseif (mean(ismember(ea_ft_senslabel('biosemi256'),         sens.label)) > 0.8)
      type = 'biosemi256';
    elseif (mean(ismember(ea_ft_senslabel('biosemi128'),         sens.label)) > 0.8)
      type = 'biosemi128';
    elseif (mean(ismember(ea_ft_senslabel('biosemi64'),          sens.label)) > 0.8)
      type = 'biosemi64';
    elseif (mean(ismember(ea_ft_senslabel('egi256'),             sens.label)) > 0.8)
      type = 'egi256';
    elseif (mean(ismember(ea_ft_senslabel('egi128'),             sens.label)) > 0.8)
      type = 'egi128';
    elseif (mean(ismember(ea_ft_senslabel('egi64'),              sens.label)) > 0.8)
      type = 'egi64';
    elseif (mean(ismember(ea_ft_senslabel('egi32'),              sens.label)) > 0.8)
      type = 'egi32';

      % the following check on the fraction of channels in the user's data rather than on the fraction of channels in the predefined set
    elseif (mean(ismember(sens.label, ea_ft_senslabel('eeg1020'))) > 0.8)
      type = 'eeg1020';
    elseif (mean(ismember(sens.label, ea_ft_senslabel('eeg1010'))) > 0.8)
      type = 'eeg1010';
    elseif (mean(ismember(sens.label, ea_ft_senslabel('eeg1005'))) > 0.8)
      type = 'eeg1005';

    elseif (sum(ismember(sens.label, ea_ft_senslabel('eeg1005'))) > 10) % Otherwise it's not even worth recognizing
      type = 'ext1020'; % this will also cover small subsets of eeg1020, eeg1010 and eeg1005
    elseif any(ismember(ea_ft_senslabel('btiref'), sens.label))
      type = 'bti'; % it might be 148 or 248 channels
    elseif any(ismember(ea_ft_senslabel('ctfref'), sens.label))
      type = 'ctf'; % it might be 151 or 275 channels

%     elseif (mean(ismember(sens.label,    ft_senslabel('nirs'))) > 0.8)
%       type = 'nirs';
    end
  end % look at label, ori and/or pnt
end % if isfield(sens, 'type')

if strcmp(type, 'unknown') && ~recursion
  % try whether only lowercase channel labels makes a difference
  if islabel
    recursion = true;
    type = ea_ft_senstype(lower(input));
    recursion = false;
  elseif isfield(input, 'label')
    input.label = lower(input.label);
    recursion = true;
    type = ea_ft_senstype(input);
    recursion = false;
  end
end

if strcmp(type, 'unknown') && ~recursion
  % try whether only uppercase channel labels makes a difference
  if islabel
    recursion = true;
    type = ea_ft_senstype(upper(input));
    recursion = false;
  elseif isfield(input, 'label')
    input.label = upper(input.label);
    recursion = true;
    type = ea_ft_senstype(input);
    recursion = false;
  end
end

if ~isempty(desired)
  % return a boolean flag
  switch desired
    case 'ext1020'
      type = any(strcmp(type, {'eeg1005' 'eeg1010' 'eeg1020' 'ext1020'}));
    case {'eeg' 'electrode'}
      type = any(strcmp(type, {'eeg' 'electrode' 'ant128' 'biosemi64' 'biosemi128' 'biosemi256' 'egi32' 'egi64' 'egi128' 'egi256' 'eeg1005' 'eeg1010' 'eeg1020' 'ext1020'}));
    case 'biosemi'
      type = any(strcmp(type, {'biosemi64' 'biosemi128' 'biosemi256'}));
    case 'egi'
      type = any(strcmp(type, {'egi32' 'egi64' 'egi128' 'egi256'}));
    case 'meg'
      type = any(strcmp(type, {'meg' 'magnetometer' 'ctf' 'bti' 'ctf64' 'ctf151' 'ctf275' 'ctf151_planar' 'ctf275_planar' 'neuromag122' 'neuromag306' 'bti148' 'bti148_planar' 'bti248' 'bti248_planar' 'bti248grad' 'bti248grad_planar' 'yokogawa9' 'yokogawa160' 'yokogawa160_planar' 'yokogawa64' 'yokogawa64_planar' 'yokogawa440' 'yokogawa440_planar' 'itab' 'itab28' 'itab153' 'itab153_planar'}));
    case 'ctf'
      type = any(strcmp(type, {'ctf' 'ctf64' 'ctf151' 'ctf275' 'ctf151_planar' 'ctf275_planar'}));
    case 'bti'
      type = any(strcmp(type, {'bti' 'bti148' 'bti148_planar' 'bti248' 'bti248_planar' 'bti248grad' 'bti248grad_planar'}));
    case 'neuromag'
      type = any(strcmp(type, {'neuromag122' 'neuromag306'}));
    case 'yokogawa'
      type = any(strcmp(type, {'yokogawa160' 'yokogawa160_planar' 'yokogawa64' 'yokogawa64_planar' 'yokogawa440' 'yokogawa440_planar'}));
    case 'itab'
      type = any(strcmp(type, {'itab' 'itab28' 'itab153' 'itab153_planar'}));
    case 'meg_axial'
      % note that neuromag306 is mixed planar and axial
      type = any(strcmp(type, {'neuromag306' 'ctf64' 'ctf151' 'ctf275' 'bti148' 'bti248' 'bti248grad' 'yokogawa9' 'yokogawa64' 'yokogawa160' 'yokogawa440'}));
    case 'meg_planar'
      % note that neuromag306 is mixed planar and axial
      type = any(strcmp(type, {'neuromag122' 'neuromag306' 'ctf151_planar' 'ctf275_planar' 'bti148_planar' 'bti248_planar' 'bti248grad_planar' 'yokogawa160_planar' 'yokogawa64_planar' 'yokogawa440_planar'}));
    otherwise
      type = any(strcmp(type, desired));
  end % switch desired
end % detemine the correspondence to the desired type

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout = {type};
previous_argin  = current_argin;
previous_argout = current_argout;

return % ft_senstype main()

function label = ea_ft_senslabel(type, varargin)

% FT_SENSLABEL returns a list of predefined sensor labels given the
% EEG or MEG system type which can be used to detect the type of data.
%
% Use as
%  label = ft_senslabel(type)
%
% The input sensor array type can be any of the following
%  'ant128'
%  'biosemi64'
%  'biosemi128'
%  'biosemi256'
%  'bti148'
%  'bti148_planar'
%  'bti248'
%  'bti248_planar'
%  'btiref'
%  'ctf151'
%  'ctf151_planar'
%  'ctf275'
%  'ctf275_planar'
%  'ctfheadloc'
%  'ctfref'
%  'eeg1005'
%  'eeg1010'
%  'eeg1020'
%  'ext1020'
%  'egi32'
%  'egi64'
%  'egi128'
%  'egi256'
%  'neuromag122'
%  'neuromag306'
%  'itab28'
%  'itab153'
%  'itab153_planar'
%  'yokogawa9'
%  'yokogawa64'
%  'yokogawa64_planar'
%  'yokogawa160'
%  'yokogawa160_planar'
%  'yokogawa440'
%  'yokogawa440_planar'
%
% It is also possible to specify
%  'eeg'
%  'electrode'
% although for these an empty set of labels (i.e. {}) will be returned.
%
% See also FT_SENSTYPE, FT_CHANNELSELECTION

% Copyright (C) 2007-2013, Robert Oostenveld
% Copyright (C) 2008, Vladimir Litvak
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%  FieldTrip is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  FieldTrip is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_senslabel.m 10340 2015-04-17 14:10:04Z jorhor $

% these are for speeding up subsequent calls with the same input arguments
persistent eeg electrode ant128 btiref bti148 bti148_planar bti148_planar_combined bti248 bti248_planar bti248_planar_combined ctfref ctfheadloc ctf64 ctf151 ctf151_planar ctf151_planar_combined ctf275 ctf275_planar ctf275_planar_combined neuromag122 neuromag122_combined neuromag306 neuromag306_combined eeg1020 eeg1010 eeg1005 ext1020 biosemi64 biosemi128 biosemi256 egi32 egi64 egi128 egi256 itab28 itab153 itab153_planar itab153_planar_combined yokogawa9 yokogawa64 yokogawa64_planar yokogawa64_planar_combined yokogawa160 yokogawa160_planar yokogawa160_planar_combined yokogawa440 yokogawa440_planar yokogawa440_planar_combined
% these are for backward compatibility
persistent neuromag122alt neuromag122alt_combined neuromag306alt neuromag306alt_combined

if nargin<1
  % ensure that all input arguments are defined
  type = 'none';
end

% get the optional input arguments
output  = ea_ft_getopt(varargin, 'output', 'normal'); % 'normal' or 'planarcombined'

if ~exist('type', 'var')
  error('the requested sensor type "%s" is not supported', type);

elseif isempty(eval(type))
  % assign the list of channels only once, keep it as persistent variable

  switch type
    case 'ant128'
      label = {
        'Z1'
        'Z2'
        'Z3'
        'Z4'
        'Z5'
        'Z6'
        'Z7'
        'Z8'
        'Z9'
        'Z10'
        'Z11'
        'Z12'
        'Z13'
        'Z14'
        'L1'
        'L2'
        'L3'
        'L4'
        'L5'
        'L6'
        'L7'
        'L8'
        'L9'
        'L10'
        'L11'
        'L12'
        'L13'
        'L14'
        'LL1'
        'LL2'
        'LL3'
        'LL4'
        'LL5'
        'LL6'
        'LL7'
        'LL8'
        'LL9'
        'LL10'
        'LL11'
        'LL12'
        'LL13'
        'LA1'
        'LA2'
        'LA3'
        'LA4'
        'LA5'
        'LB1'
        'LB2'
        'LB3'
        'LB4'
        'LB5'
        'LB6'
        'LC1'
        'LC2'
        'LC3'
        'LC4'
        'LC5'
        'LC6'
        'LC7'
        'LD1'
        'LD2'
        'LD3'
        'LD4'
        'LD5'
        'LD6'
        'LD7'
        'LE1'
        'LE2'
        'LE3'
        'LE4'
        'Lm'
        'R1'
        'R2'
        'R3'
        'R4'
        'R5'
        'R6'
        'R7'
        'R8'
        'R9'
        'R10'
        'R11'
        'R12'
        'R13'
        'R14'
        'RR1'
        'RR2'
        'RR3'
        'RR4'
        'RR5'
        'RR6'
        'RR7'
        'RR8'
        'RR9'
        'RR10'
        'RR11'
        'RR12'
        'RR13'
        'RA1'
        'RA2'
        'RA3'
        'RA4'
        'RA5'
        'RB1'
        'RB2'
        'RB3'
        'RB4'
        'RB5'
        'RB6'
        'RC1'
        'RC2'
        'RC3'
        'RC4'
        'RC5'
        'RC6'
        'RC7'
        'RD1'
        'RD2'
        'RD3'
        'RD4'
        'RD5'
        'RD6'
        'RD7'
        'RE1'
        'RE2'
        'RE3'
        'RE4'
        'Rm'
        };

    case 'btiref'
      label = {
        'MRxA'
        'MRyA'
        'MRzA'
        'MLxA'
        'MLyA'
        'MLzA'
        'MCxA'
        'MCyA'
        'MCzA'
        'MRxaA'
        'MRyaA'
        'MRzaA'
        'MLxaA'
        'MLyaA'
        'MLzaA'
        'MCxaA'
        'MCyaA'
        'MCzaA'
        'GxxA'
        'GyxA'
        'GzxA'
        'GyyA'
        'GzyA'
        };

    case 'bti148'
      label = cell(148,1);
      for i=1:148
        label{i,1} = sprintf('A%d', i);
      end

    case 'bti148_planar'
      label = cell(148,3);
      for i=1:148
        label{i,1} = sprintf('A%d_dH', i);
        label{i,2} = sprintf('A%d_dV', i);
        label{i,3} = sprintf('A%d', i);
      end
      bti148_planar_combined = label(:,3);
      label = label(:,1:2);

    case 'bti248'
      label = cell(248,1);
      for i=1:248
        label{i,1} = sprintf('A%d', i);
      end

    case 'bti248_planar'
      label = cell(248,3);
      for i=1:248
        label{i,1} = sprintf('A%d_dH', i);
        label{i,2} = sprintf('A%d_dV', i);
        label{i,3} = sprintf('A%d', i);
      end
      bti248_planar_combined = label(:,3);
      label = label(:,1:2);

    case 'ctfref'
      label = {
        'BG1'
        'BG2'
        'BG3'
        'BP1'
        'BP2'
        'BP3'
        'BR1'
        'BR2'
        'BR3'
        'G11'
        'G12'
        'G13'
        'G22'
        'G23'
        'P11'
        'P12'
        'P13'
        'P22'
        'P23'
        'Q11'
        'Q12'
        'Q13'
        'Q22'
        'Q23'
        'R11'
        'R12'
        'R13'
        'R22'
        'R23'
        };

    case 'ctfheadloc'
      label = {
        'HLC0011'
        'HLC0012'
        'HLC0013'
        'HLC0021'
        'HLC0022'
        'HLC0023'
        'HLC0031'
        'HLC0032'
        'HLC0033'
        'HLC0018'
        'HLC0028'
        'HLC0038'
        'HLC0014'
        'HLC0015'
        'HLC0016'
        'HLC0017'
        'HLC0024'
        'HLC0025'
        'HLC0026'
        'HLC0027'
        'HLC0034'
        'HLC0035'
        'HLC0036'
        'HLC0037'
        };

    case 'ctf64'
      label = {
        'SL11'
        'SL12'
        'SL13'
        'SL14'
        'SL15'
        'SL16'
        'SL17'
        'SL18'
        'SL19'
        'SL21'
        'SL22'
        'SL23'
        'SL24'
        'SL25'
        'SL26'
        'SL27'
        'SL28'
        'SL29'
        'SL31'
        'SL32'
        'SL33'
        'SL34'
        'SL35'
        'SL41'
        'SL42'
        'SL43'
        'SL44'
        'SL45'
        'SL46'
        'SL47'
        'SL51'
        'SL52'
        'SR11'
        'SR12'
        'SR13'
        'SR14'
        'SR15'
        'SR16'
        'SR17'
        'SR18'
        'SR19'
        'SR21'
        'SR22'
        'SR23'
        'SR24'
        'SR25'
        'SR26'
        'SR27'
        'SR28'
        'SR29'
        'SR31'
        'SR32'
        'SR33'
        'SR34'
        'SR35'
        'SR41'
        'SR42'
        'SR43'
        'SR44'
        'SR45'
        'SR46'
        'SR47'
        'SR51'
        'SR52'
        };

    case 'ctf151'
      label = {
        'MLC11'
        'MLC12'
        'MLC13'
        'MLC14'
        'MLC15'
        'MLC21'
        'MLC22'
        'MLC23'
        'MLC24'
        'MLC31'
        'MLC32'
        'MLC33'
        'MLC41'
        'MLC42'
        'MLC43'
        'MLF11'
        'MLF12'
        'MLF21'
        'MLF22'
        'MLF23'
        'MLF31'
        'MLF32'
        'MLF33'
        'MLF34'
        'MLF41'
        'MLF42'
        'MLF43'
        'MLF44'
        'MLF45'
        'MLF51'
        'MLF52'
        'MLO11'
        'MLO12'
        'MLO21'
        'MLO22'
        'MLO31'
        'MLO32'
        'MLO33'
        'MLO41'
        'MLO42'
        'MLO43'
        'MLP11'
        'MLP12'
        'MLP13'
        'MLP21'
        'MLP22'
        'MLP31'
        'MLP32'
        'MLP33'
        'MLP34'
        'MLT11'
        'MLT12'
        'MLT13'
        'MLT14'
        'MLT15'
        'MLT16'
        'MLT21'
        'MLT22'
        'MLT23'
        'MLT24'
        'MLT25'
        'MLT26'
        'MLT31'
        'MLT32'
        'MLT33'
        'MLT34'
        'MLT35'
        'MLT41'
        'MLT42'
        'MLT43'
        'MLT44'
        'MRC11'
        'MRC12'
        'MRC13'
        'MRC14'
        'MRC15'
        'MRC21'
        'MRC22'
        'MRC23'
        'MRC24'
        'MRC31'
        'MRC32'
        'MRC33'
        'MRC41'
        'MRC42'
        'MRC43'
        'MRF11'
        'MRF12'
        'MRF21'
        'MRF22'
        'MRF23'
        'MRF31'
        'MRF32'
        'MRF33'
        'MRF34'
        'MRF41'
        'MRF42'
        'MRF43'
        'MRF44'
        'MRF45'
        'MRF51'
        'MRF52'
        'MRO11'
        'MRO12'
        'MRO21'
        'MRO22'
        'MRO31'
        'MRO32'
        'MRO33'
        'MRO41'
        'MRO42'
        'MRO43'
        'MRP11'
        'MRP12'
        'MRP13'
        'MRP21'
        'MRP22'
        'MRP31'
        'MRP32'
        'MRP33'
        'MRP34'
        'MRT11'
        'MRT12'
        'MRT13'
        'MRT14'
        'MRT15'
        'MRT16'
        'MRT21'
        'MRT22'
        'MRT23'
        'MRT24'
        'MRT25'
        'MRT26'
        'MRT31'
        'MRT32'
        'MRT33'
        'MRT34'
        'MRT35'
        'MRT41'
        'MRT42'
        'MRT43'
        'MRT44'
        'MZC01'
        'MZC02'
        'MZF01'
        'MZF02'
        'MZF03'
        'MZO01'
        'MZO02'
        'MZP01'
        'MZP02'
        };

    case 'ctf151_planar'
      label = {
        'MLC11_dH'  'MLC11_dV'  'MLC11'
        'MLC12_dH'  'MLC12_dV'  'MLC12'
        'MLC13_dH'  'MLC13_dV'  'MLC13'
        'MLC14_dH'  'MLC14_dV'  'MLC14'
        'MLC15_dH'  'MLC15_dV'  'MLC15'
        'MLC21_dH'  'MLC21_dV'  'MLC21'
        'MLC22_dH'  'MLC22_dV'  'MLC22'
        'MLC23_dH'  'MLC23_dV'  'MLC23'
        'MLC24_dH'  'MLC24_dV'  'MLC24'
        'MLC31_dH'  'MLC31_dV'  'MLC31'
        'MLC32_dH'  'MLC32_dV'  'MLC32'
        'MLC33_dH'  'MLC33_dV'  'MLC33'
        'MLC41_dH'  'MLC41_dV'  'MLC41'
        'MLC42_dH'  'MLC42_dV'  'MLC42'
        'MLC43_dH'  'MLC43_dV'  'MLC43'
        'MLF11_dH'  'MLF11_dV'  'MLF11'
        'MLF12_dH'  'MLF12_dV'  'MLF12'
        'MLF21_dH'  'MLF21_dV'  'MLF21'
        'MLF22_dH'  'MLF22_dV'  'MLF22'
        'MLF23_dH'  'MLF23_dV'  'MLF23'
        'MLF31_dH'  'MLF31_dV'  'MLF31'
        'MLF32_dH'  'MLF32_dV'  'MLF32'
        'MLF33_dH'  'MLF33_dV'  'MLF33'
        'MLF34_dH'  'MLF34_dV'  'MLF34'
        'MLF41_dH'  'MLF41_dV'  'MLF41'
        'MLF42_dH'  'MLF42_dV'  'MLF42'
        'MLF43_dH'  'MLF43_dV'  'MLF43'
        'MLF44_dH'  'MLF44_dV'  'MLF44'
        'MLF45_dH'  'MLF45_dV'  'MLF45'
        'MLF51_dH'  'MLF51_dV'  'MLF51'
        'MLF52_dH'  'MLF52_dV'  'MLF52'
        'MLO11_dH'  'MLO11_dV'  'MLO11'
        'MLO12_dH'  'MLO12_dV'  'MLO12'
        'MLO21_dH'  'MLO21_dV'  'MLO21'
        'MLO22_dH'  'MLO22_dV'  'MLO22'
        'MLO31_dH'  'MLO31_dV'  'MLO31'
        'MLO32_dH'  'MLO32_dV'  'MLO32'
        'MLO33_dH'  'MLO33_dV'  'MLO33'
        'MLO41_dH'  'MLO41_dV'  'MLO41'
        'MLO42_dH'  'MLO42_dV'  'MLO42'
        'MLO43_dH'  'MLO43_dV'  'MLO43'
        'MLP11_dH'  'MLP11_dV'  'MLP11'
        'MLP12_dH'  'MLP12_dV'  'MLP12'
        'MLP13_dH'  'MLP13_dV'  'MLP13'
        'MLP21_dH'  'MLP21_dV'  'MLP21'
        'MLP22_dH'  'MLP22_dV'  'MLP22'
        'MLP31_dH'  'MLP31_dV'  'MLP31'
        'MLP32_dH'  'MLP32_dV'  'MLP32'
        'MLP33_dH'  'MLP33_dV'  'MLP33'
        'MLP34_dH'  'MLP34_dV'  'MLP34'
        'MLT11_dH'  'MLT11_dV'  'MLT11'
        'MLT12_dH'  'MLT12_dV'  'MLT12'
        'MLT13_dH'  'MLT13_dV'  'MLT13'
        'MLT14_dH'  'MLT14_dV'  'MLT14'
        'MLT15_dH'  'MLT15_dV'  'MLT15'
        'MLT16_dH'  'MLT16_dV'  'MLT16'
        'MLT21_dH'  'MLT21_dV'  'MLT21'
        'MLT22_dH'  'MLT22_dV'  'MLT22'
        'MLT23_dH'  'MLT23_dV'  'MLT23'
        'MLT24_dH'  'MLT24_dV'  'MLT24'
        'MLT25_dH'  'MLT25_dV'  'MLT25'
        'MLT26_dH'  'MLT26_dV'  'MLT26'
        'MLT31_dH'  'MLT31_dV'  'MLT31'
        'MLT32_dH'  'MLT32_dV'  'MLT32'
        'MLT33_dH'  'MLT33_dV'  'MLT33'
        'MLT34_dH'  'MLT34_dV'  'MLT34'
        'MLT35_dH'  'MLT35_dV'  'MLT35'
        'MLT41_dH'  'MLT41_dV'  'MLT41'
        'MLT42_dH'  'MLT42_dV'  'MLT42'
        'MLT43_dH'  'MLT43_dV'  'MLT43'
        'MLT44_dH'  'MLT44_dV'  'MLT44'
        'MRC11_dH'  'MRC11_dV'  'MRC11'
        'MRC12_dH'  'MRC12_dV'  'MRC12'
        'MRC13_dH'  'MRC13_dV'  'MRC13'
        'MRC14_dH'  'MRC14_dV'  'MRC14'
        'MRC15_dH'  'MRC15_dV'  'MRC15'
        'MRC21_dH'  'MRC21_dV'  'MRC21'
        'MRC22_dH'  'MRC22_dV'  'MRC22'
        'MRC23_dH'  'MRC23_dV'  'MRC23'
        'MRC24_dH'  'MRC24_dV'  'MRC24'
        'MRC31_dH'  'MRC31_dV'  'MRC31'
        'MRC32_dH'  'MRC32_dV'  'MRC32'
        'MRC33_dH'  'MRC33_dV'  'MRC33'
        'MRC41_dH'  'MRC41_dV'  'MRC41'
        'MRC42_dH'  'MRC42_dV'  'MRC42'
        'MRC43_dH'  'MRC43_dV'  'MRC43'
        'MRF11_dH'  'MRF11_dV'  'MRF11'
        'MRF12_dH'  'MRF12_dV'  'MRF12'
        'MRF21_dH'  'MRF21_dV'  'MRF21'
        'MRF22_dH'  'MRF22_dV'  'MRF22'
        'MRF23_dH'  'MRF23_dV'  'MRF23'
        'MRF31_dH'  'MRF31_dV'  'MRF31'
        'MRF32_dH'  'MRF32_dV'  'MRF32'
        'MRF33_dH'  'MRF33_dV'  'MRF33'
        'MRF34_dH'  'MRF34_dV'  'MRF34'
        'MRF41_dH'  'MRF41_dV'  'MRF41'
        'MRF42_dH'  'MRF42_dV'  'MRF42'
        'MRF43_dH'  'MRF43_dV'  'MRF43'
        'MRF44_dH'  'MRF44_dV'  'MRF44'
        'MRF45_dH'  'MRF45_dV'  'MRF45'
        'MRF51_dH'  'MRF51_dV'  'MRF51'
        'MRF52_dH'  'MRF52_dV'  'MRF52'
        'MRO11_dH'  'MRO11_dV'  'MRO11'
        'MRO12_dH'  'MRO12_dV'  'MRO12'
        'MRO21_dH'  'MRO21_dV'  'MRO21'
        'MRO22_dH'  'MRO22_dV'  'MRO22'
        'MRO31_dH'  'MRO31_dV'  'MRO31'
        'MRO32_dH'  'MRO32_dV'  'MRO32'
        'MRO33_dH'  'MRO33_dV'  'MRO33'
        'MRO41_dH'  'MRO41_dV'  'MRO41'
        'MRO42_dH'  'MRO42_dV'  'MRO42'
        'MRO43_dH'  'MRO43_dV'  'MRO43'
        'MRP11_dH'  'MRP11_dV'  'MRP11'
        'MRP12_dH'  'MRP12_dV'  'MRP12'
        'MRP13_dH'  'MRP13_dV'  'MRP13'
        'MRP21_dH'  'MRP21_dV'  'MRP21'
        'MRP22_dH'  'MRP22_dV'  'MRP22'
        'MRP31_dH'  'MRP31_dV'  'MRP31'
        'MRP32_dH'  'MRP32_dV'  'MRP32'
        'MRP33_dH'  'MRP33_dV'  'MRP33'
        'MRP34_dH'  'MRP34_dV'  'MRP34'
        'MRT11_dH'  'MRT11_dV'  'MRT11'
        'MRT12_dH'  'MRT12_dV'  'MRT12'
        'MRT13_dH'  'MRT13_dV'  'MRT13'
        'MRT14_dH'  'MRT14_dV'  'MRT14'
        'MRT15_dH'  'MRT15_dV'  'MRT15'
        'MRT16_dH'  'MRT16_dV'  'MRT16'
        'MRT21_dH'  'MRT21_dV'  'MRT21'
        'MRT22_dH'  'MRT22_dV'  'MRT22'
        'MRT23_dH'  'MRT23_dV'  'MRT23'
        'MRT24_dH'  'MRT24_dV'  'MRT24'
        'MRT25_dH'  'MRT25_dV'  'MRT25'
        'MRT26_dH'  'MRT26_dV'  'MRT26'
        'MRT31_dH'  'MRT31_dV'  'MRT31'
        'MRT32_dH'  'MRT32_dV'  'MRT32'
        'MRT33_dH'  'MRT33_dV'  'MRT33'
        'MRT34_dH'  'MRT34_dV'  'MRT34'
        'MRT35_dH'  'MRT35_dV'  'MRT35'
        'MRT41_dH'  'MRT41_dV'  'MRT41'
        'MRT42_dH'  'MRT42_dV'  'MRT42'
        'MRT43_dH'  'MRT43_dV'  'MRT43'
        'MRT44_dH'  'MRT44_dV'  'MRT44'
        'MZC01_dH'  'MZC01_dV'  'MZC01'
        'MZC02_dH'  'MZC02_dV'  'MZC02'
        'MZF01_dH'  'MZF01_dV'  'MZF01'
        'MZF02_dH'  'MZF02_dV'  'MZF02'
        'MZF03_dH'  'MZF03_dV'  'MZF03'
        'MZO01_dH'  'MZO01_dV'  'MZO01'
        'MZO02_dH'  'MZO02_dV'  'MZO02'
        'MZP01_dH'  'MZP01_dV'  'MZP01'
        'MZP02_dH'  'MZP02_dV'  'MZP02'
        };
      ctf151_planar_combined = label(:,3);
      label = label(:,1:2);

    case 'ctf275'
      label = {
        'MLC11'
        'MLC12'
        'MLC13'
        'MLC14'
        'MLC15'
        'MLC16'
        'MLC17'
        'MLC21'
        'MLC22'
        'MLC23'
        'MLC24'
        'MLC25'
        'MLC31'
        'MLC32'
        'MLC41'
        'MLC42'
        'MLC51'
        'MLC52'
        'MLC53'
        'MLC54'
        'MLC55'
        'MLC61'
        'MLC62'
        'MLC63'
        'MLF11'
        'MLF12'
        'MLF13'
        'MLF14'
        'MLF21'
        'MLF22'
        'MLF23'
        'MLF24'
        'MLF25'
        'MLF31'
        'MLF32'
        'MLF33'
        'MLF34'
        'MLF35'
        'MLF41'
        'MLF42'
        'MLF43'
        'MLF44'
        'MLF45'
        'MLF46'
        'MLF51'
        'MLF52'
        'MLF53'
        'MLF54'
        'MLF55'
        'MLF56'
        'MLF61'
        'MLF62'
        'MLF63'
        'MLF64'
        'MLF65'
        'MLF66'
        'MLF67'
        'MLO11'
        'MLO12'
        'MLO13'
        'MLO14'
        'MLO21'
        'MLO22'
        'MLO23'
        'MLO24'
        'MLO31'
        'MLO32'
        'MLO33'
        'MLO34'
        'MLO41'
        'MLO42'
        'MLO43'
        'MLO44'
        'MLO51'
        'MLO52'
        'MLO53'
        'MLP11'
        'MLP12'
        'MLP21'
        'MLP22'
        'MLP23'
        'MLP31'
        'MLP32'
        'MLP33'
        'MLP34'
        'MLP35'
        'MLP41'
        'MLP42'
        'MLP43'
        'MLP44'
        'MLP45'
        'MLP51'
        'MLP52'
        'MLP53'
        'MLP54'
        'MLP55'
        'MLP56'
        'MLP57'
        'MLT11'
        'MLT12'
        'MLT13'
        'MLT14'
        'MLT15'
        'MLT16'
        'MLT21'
        'MLT22'
        'MLT23'
        'MLT24'
        'MLT25'
        'MLT26'
        'MLT27'
        'MLT31'
        'MLT32'
        'MLT33'
        'MLT34'
        'MLT35'
        'MLT36'
        'MLT37'
        'MLT41'
        'MLT42'
        'MLT43'
        'MLT44'
        'MLT45'
        'MLT46'
        'MLT47'
        'MLT51'
        'MLT52'
        'MLT53'
        'MLT54'
        'MLT55'
        'MLT56'
        'MLT57'
        'MRC11'
        'MRC12'
        'MRC13'
        'MRC14'
        'MRC15'
        'MRC16'
        'MRC17'
        'MRC21'
        'MRC22'
        'MRC23'
        'MRC24'
        'MRC25'
        'MRC31'
        'MRC32'
        'MRC41'
        'MRC42'
        'MRC51'
        'MRC52'
        'MRC53'
        'MRC54'
        'MRC55'
        'MRC61'
        'MRC62'
        'MRC63'
        'MRF11'
        'MRF12'
        'MRF13'
        'MRF14'
        'MRF21'
        'MRF22'
        'MRF23'
        'MRF24'
        'MRF25'
        'MRF31'
        'MRF32'
        'MRF33'
        'MRF34'
        'MRF35'
        'MRF41'
        'MRF42'
        'MRF43'
        'MRF44'
        'MRF45'
        'MRF46'
        'MRF51'
        'MRF52'
        'MRF53'
        'MRF54'
        'MRF55'
        'MRF56'
        'MRF61'
        'MRF62'
        'MRF63'
        'MRF64'
        'MRF65'
        'MRF66'
        'MRF67'
        'MRO11'
        'MRO12'
        'MRO13'
        'MRO14'
        'MRO21'
        'MRO22'
        'MRO23'
        'MRO24'
        'MRO31'
        'MRO32'
        'MRO33'
        'MRO34'
        'MRO41'
        'MRO42'
        'MRO43'
        'MRO44'
        'MRO51'
        'MRO52'
        'MRO53'
        'MRP11'
        'MRP12'
        'MRP21'
        'MRP22'
        'MRP23'
        'MRP31'
        'MRP32'
        'MRP33'
        'MRP34'
        'MRP35'
        'MRP41'
        'MRP42'
        'MRP43'
        'MRP44'
        'MRP45'
        'MRP51'
        'MRP52'
        'MRP53'
        'MRP54'
        'MRP55'
        'MRP56'
        'MRP57'
        'MRT11'
        'MRT12'
        'MRT13'
        'MRT14'
        'MRT15'
        'MRT16'
        'MRT21'
        'MRT22'
        'MRT23'
        'MRT24'
        'MRT25'
        'MRT26'
        'MRT27'
        'MRT31'
        'MRT32'
        'MRT33'
        'MRT34'
        'MRT35'
        'MRT36'
        'MRT37'
        'MRT41'
        'MRT42'
        'MRT43'
        'MRT44'
        'MRT45'
        'MRT46'
        'MRT47'
        'MRT51'
        'MRT52'
        'MRT53'
        'MRT54'
        'MRT55'
        'MRT56'
        'MRT57'
        'MZC01'
        'MZC02'
        'MZC03'
        'MZC04'
        'MZF01'
        'MZF02'
        'MZF03'
        'MZO01'
        'MZO02'
        'MZO03'
        'MZP01'
        };

    case 'ctf275_planar'
      label = {
        'MLC11_dH'  'MLC11_dV'  'MLC11'
        'MLC12_dH'  'MLC12_dV'  'MLC12'
        'MLC13_dH'  'MLC13_dV'  'MLC13'
        'MLC14_dH'  'MLC14_dV'  'MLC14'
        'MLC15_dH'  'MLC15_dV'  'MLC15'
        'MLC16_dH'  'MLC16_dV'  'MLC16'
        'MLC17_dH'  'MLC17_dV'  'MLC17'
        'MLC21_dH'  'MLC21_dV'  'MLC21'
        'MLC22_dH'  'MLC22_dV'  'MLC22'
        'MLC23_dH'  'MLC23_dV'  'MLC23'
        'MLC24_dH'  'MLC24_dV'  'MLC24'
        'MLC25_dH'  'MLC25_dV'  'MLC25'
        'MLC31_dH'  'MLC31_dV'  'MLC31'
        'MLC32_dH'  'MLC32_dV'  'MLC32'
        'MLC41_dH'  'MLC41_dV'  'MLC41'
        'MLC42_dH'  'MLC42_dV'  'MLC42'
        'MLC51_dH'  'MLC51_dV'  'MLC51'
        'MLC52_dH'  'MLC52_dV'  'MLC52'
        'MLC53_dH'  'MLC53_dV'  'MLC53'
        'MLC54_dH'  'MLC54_dV'  'MLC54'
        'MLC55_dH'  'MLC55_dV'  'MLC55'
        'MLC61_dH'  'MLC61_dV'  'MLC61'
        'MLC62_dH'  'MLC62_dV'  'MLC62'
        'MLC63_dH'  'MLC63_dV'  'MLC63'
        'MLF11_dH'  'MLF11_dV'  'MLF11'
        'MLF12_dH'  'MLF12_dV'  'MLF12'
        'MLF13_dH'  'MLF13_dV'  'MLF13'
        'MLF14_dH'  'MLF14_dV'  'MLF14'
        'MLF21_dH'  'MLF21_dV'  'MLF21'
        'MLF22_dH'  'MLF22_dV'  'MLF22'
        'MLF23_dH'  'MLF23_dV'  'MLF23'
        'MLF24_dH'  'MLF24_dV'  'MLF24'
        'MLF25_dH'  'MLF25_dV'  'MLF25'
        'MLF31_dH'  'MLF31_dV'  'MLF31'
        'MLF32_dH'  'MLF32_dV'  'MLF32'
        'MLF33_dH'  'MLF33_dV'  'MLF33'
        'MLF34_dH'  'MLF34_dV'  'MLF34'
        'MLF35_dH'  'MLF35_dV'  'MLF35'
        'MLF41_dH'  'MLF41_dV'  'MLF41'
        'MLF42_dH'  'MLF42_dV'  'MLF42'
        'MLF43_dH'  'MLF43_dV'  'MLF43'
        'MLF44_dH'  'MLF44_dV'  'MLF44'
        'MLF45_dH'  'MLF45_dV'  'MLF45'
        'MLF46_dH'  'MLF46_dV'  'MLF46'
        'MLF51_dH'  'MLF51_dV'  'MLF51'
        'MLF52_dH'  'MLF52_dV'  'MLF52'
        'MLF53_dH'  'MLF53_dV'  'MLF53'
        'MLF54_dH'  'MLF54_dV'  'MLF54'
        'MLF55_dH'  'MLF55_dV'  'MLF55'
        'MLF56_dH'  'MLF56_dV'  'MLF56'
        'MLF61_dH'  'MLF61_dV'  'MLF61'
        'MLF62_dH'  'MLF62_dV'  'MLF62'
        'MLF63_dH'  'MLF63_dV'  'MLF63'
        'MLF64_dH'  'MLF64_dV'  'MLF64'
        'MLF65_dH'  'MLF65_dV'  'MLF65'
        'MLF66_dH'  'MLF66_dV'  'MLF66'
        'MLF67_dH'  'MLF67_dV'  'MLF67'
        'MLO11_dH'  'MLO11_dV'  'MLO11'
        'MLO12_dH'  'MLO12_dV'  'MLO12'
        'MLO13_dH'  'MLO13_dV'  'MLO13'
        'MLO14_dH'  'MLO14_dV'  'MLO14'
        'MLO21_dH'  'MLO21_dV'  'MLO21'
        'MLO22_dH'  'MLO22_dV'  'MLO22'
        'MLO23_dH'  'MLO23_dV'  'MLO23'
        'MLO24_dH'  'MLO24_dV'  'MLO24'
        'MLO31_dH'  'MLO31_dV'  'MLO31'
        'MLO32_dH'  'MLO32_dV'  'MLO32'
        'MLO33_dH'  'MLO33_dV'  'MLO33'
        'MLO34_dH'  'MLO34_dV'  'MLO34'
        'MLO41_dH'  'MLO41_dV'  'MLO41'
        'MLO42_dH'  'MLO42_dV'  'MLO42'
        'MLO43_dH'  'MLO43_dV'  'MLO43'
        'MLO44_dH'  'MLO44_dV'  'MLO44'
        'MLO51_dH'  'MLO51_dV'  'MLO51'
        'MLO52_dH'  'MLO52_dV'  'MLO52'
        'MLO53_dH'  'MLO53_dV'  'MLO53'
        'MLP11_dH'  'MLP11_dV'  'MLP11'
        'MLP12_dH'  'MLP12_dV'  'MLP12'
        'MLP21_dH'  'MLP21_dV'  'MLP21'
        'MLP22_dH'  'MLP22_dV'  'MLP22'
        'MLP23_dH'  'MLP23_dV'  'MLP23'
        'MLP31_dH'  'MLP31_dV'  'MLP31'
        'MLP32_dH'  'MLP32_dV'  'MLP32'
        'MLP33_dH'  'MLP33_dV'  'MLP33'
        'MLP34_dH'  'MLP34_dV'  'MLP34'
        'MLP35_dH'  'MLP35_dV'  'MLP35'
        'MLP41_dH'  'MLP41_dV'  'MLP41'
        'MLP42_dH'  'MLP42_dV'  'MLP42'
        'MLP43_dH'  'MLP43_dV'  'MLP43'
        'MLP44_dH'  'MLP44_dV'  'MLP44'
        'MLP45_dH'  'MLP45_dV'  'MLP45'
        'MLP51_dH'  'MLP51_dV'  'MLP51'
        'MLP52_dH'  'MLP52_dV'  'MLP52'
        'MLP53_dH'  'MLP53_dV'  'MLP53'
        'MLP54_dH'  'MLP54_dV'  'MLP54'
        'MLP55_dH'  'MLP55_dV'  'MLP55'
        'MLP56_dH'  'MLP56_dV'  'MLP56'
        'MLP57_dH'  'MLP57_dV'  'MLP57'
        'MLT11_dH'  'MLT11_dV'  'MLT11'
        'MLT12_dH'  'MLT12_dV'  'MLT12'
        'MLT13_dH'  'MLT13_dV'  'MLT13'
        'MLT14_dH'  'MLT14_dV'  'MLT14'
        'MLT15_dH'  'MLT15_dV'  'MLT15'
        'MLT16_dH'  'MLT16_dV'  'MLT16'
        'MLT21_dH'  'MLT21_dV'  'MLT21'
        'MLT22_dH'  'MLT22_dV'  'MLT22'
        'MLT23_dH'  'MLT23_dV'  'MLT23'
        'MLT24_dH'  'MLT24_dV'  'MLT24'
        'MLT25_dH'  'MLT25_dV'  'MLT25'
        'MLT26_dH'  'MLT26_dV'  'MLT26'
        'MLT27_dH'  'MLT27_dV'  'MLT27'
        'MLT31_dH'  'MLT31_dV'  'MLT31'
        'MLT32_dH'  'MLT32_dV'  'MLT32'
        'MLT33_dH'  'MLT33_dV'  'MLT33'
        'MLT34_dH'  'MLT34_dV'  'MLT34'
        'MLT35_dH'  'MLT35_dV'  'MLT35'
        'MLT36_dH'  'MLT36_dV'  'MLT36'
        'MLT37_dH'  'MLT37_dV'  'MLT37'
        'MLT41_dH'  'MLT41_dV'  'MLT41'
        'MLT42_dH'  'MLT42_dV'  'MLT42'
        'MLT43_dH'  'MLT43_dV'  'MLT43'
        'MLT44_dH'  'MLT44_dV'  'MLT44'
        'MLT45_dH'  'MLT45_dV'  'MLT45'
        'MLT46_dH'  'MLT46_dV'  'MLT46'
        'MLT47_dH'  'MLT47_dV'  'MLT47'
        'MLT51_dH'  'MLT51_dV'  'MLT51'
        'MLT52_dH'  'MLT52_dV'  'MLT52'
        'MLT53_dH'  'MLT53_dV'  'MLT53'
        'MLT54_dH'  'MLT54_dV'  'MLT54'
        'MLT55_dH'  'MLT55_dV'  'MLT55'
        'MLT56_dH'  'MLT56_dV'  'MLT56'
        'MLT57_dH'  'MLT57_dV'  'MLT57'
        'MRC11_dH'  'MRC11_dV'  'MRC11'
        'MRC12_dH'  'MRC12_dV'  'MRC12'
        'MRC13_dH'  'MRC13_dV'  'MRC13'
        'MRC14_dH'  'MRC14_dV'  'MRC14'
        'MRC15_dH'  'MRC15_dV'  'MRC15'
        'MRC16_dH'  'MRC16_dV'  'MRC16'
        'MRC17_dH'  'MRC17_dV'  'MRC17'
        'MRC21_dH'  'MRC21_dV'  'MRC21'
        'MRC22_dH'  'MRC22_dV'  'MRC22'
        'MRC23_dH'  'MRC23_dV'  'MRC23'
        'MRC24_dH'  'MRC24_dV'  'MRC24'
        'MRC25_dH'  'MRC25_dV'  'MRC25'
        'MRC31_dH'  'MRC31_dV'  'MRC31'
        'MRC32_dH'  'MRC32_dV'  'MRC32'
        'MRC41_dH'  'MRC41_dV'  'MRC41'
        'MRC42_dH'  'MRC42_dV'  'MRC42'
        'MRC51_dH'  'MRC51_dV'  'MRC51'
        'MRC52_dH'  'MRC52_dV'  'MRC52'
        'MRC53_dH'  'MRC53_dV'  'MRC53'
        'MRC54_dH'  'MRC54_dV'  'MRC54'
        'MRC55_dH'  'MRC55_dV'  'MRC55'
        'MRC61_dH'  'MRC61_dV'  'MRC61'
        'MRC62_dH'  'MRC62_dV'  'MRC62'
        'MRC63_dH'  'MRC63_dV'  'MRC63'
        'MRF11_dH'  'MRF11_dV'  'MRF11'
        'MRF12_dH'  'MRF12_dV'  'MRF12'
        'MRF13_dH'  'MRF13_dV'  'MRF13'
        'MRF14_dH'  'MRF14_dV'  'MRF14'
        'MRF21_dH'  'MRF21_dV'  'MRF21'
        'MRF22_dH'  'MRF22_dV'  'MRF22'
        'MRF23_dH'  'MRF23_dV'  'MRF23'
        'MRF24_dH'  'MRF24_dV'  'MRF24'
        'MRF25_dH'  'MRF25_dV'  'MRF25'
        'MRF31_dH'  'MRF31_dV'  'MRF31'
        'MRF32_dH'  'MRF32_dV'  'MRF32'
        'MRF33_dH'  'MRF33_dV'  'MRF33'
        'MRF34_dH'  'MRF34_dV'  'MRF34'
        'MRF35_dH'  'MRF35_dV'  'MRF35'
        'MRF41_dH'  'MRF41_dV'  'MRF41'
        'MRF42_dH'  'MRF42_dV'  'MRF42'
        'MRF43_dH'  'MRF43_dV'  'MRF43'
        'MRF44_dH'  'MRF44_dV'  'MRF44'
        'MRF45_dH'  'MRF45_dV'  'MRF45'
        'MRF46_dH'  'MRF46_dV'  'MRF46'
        'MRF51_dH'  'MRF51_dV'  'MRF51'
        'MRF52_dH'  'MRF52_dV'  'MRF52'
        'MRF53_dH'  'MRF53_dV'  'MRF53'
        'MRF54_dH'  'MRF54_dV'  'MRF54'
        'MRF55_dH'  'MRF55_dV'  'MRF55'
        'MRF56_dH'  'MRF56_dV'  'MRF56'
        'MRF61_dH'  'MRF61_dV'  'MRF61'
        'MRF62_dH'  'MRF62_dV'  'MRF62'
        'MRF63_dH'  'MRF63_dV'  'MRF63'
        'MRF64_dH'  'MRF64_dV'  'MRF64'
        'MRF65_dH'  'MRF65_dV'  'MRF65'
        'MRF66_dH'  'MRF66_dV'  'MRF66'
        'MRF67_dH'  'MRF67_dV'  'MRF67'
        'MRO11_dH'  'MRO11_dV'  'MRO11'
        'MRO12_dH'  'MRO12_dV'  'MRO12'
        'MRO13_dH'  'MRO13_dV'  'MRO13'
        'MRO14_dH'  'MRO14_dV'  'MRO14'
        'MRO21_dH'  'MRO21_dV'  'MRO21'
        'MRO22_dH'  'MRO22_dV'  'MRO22'
        'MRO23_dH'  'MRO23_dV'  'MRO23'
        'MRO24_dH'  'MRO24_dV'  'MRO24'
        'MRO31_dH'  'MRO31_dV'  'MRO31'
        'MRO32_dH'  'MRO32_dV'  'MRO32'
        'MRO33_dH'  'MRO33_dV'  'MRO33'
        'MRO34_dH'  'MRO34_dV'  'MRO34'
        'MRO41_dH'  'MRO41_dV'  'MRO41'
        'MRO42_dH'  'MRO42_dV'  'MRO42'
        'MRO43_dH'  'MRO43_dV'  'MRO43'
        'MRO44_dH'  'MRO44_dV'  'MRO44'
        'MRO51_dH'  'MRO51_dV'  'MRO51'
        'MRO52_dH'  'MRO52_dV'  'MRO52'
        'MRO53_dH'  'MRO53_dV'  'MRO53'
        'MRP11_dH'  'MRP11_dV'  'MRP11'
        'MRP12_dH'  'MRP12_dV'  'MRP12'
        'MRP21_dH'  'MRP21_dV'  'MRP21'
        'MRP22_dH'  'MRP22_dV'  'MRP22'
        'MRP23_dH'  'MRP23_dV'  'MRP23'
        'MRP31_dH'  'MRP31_dV'  'MRP31'
        'MRP32_dH'  'MRP32_dV'  'MRP32'
        'MRP33_dH'  'MRP33_dV'  'MRP33'
        'MRP34_dH'  'MRP34_dV'  'MRP34'
        'MRP35_dH'  'MRP35_dV'  'MRP35'
        'MRP41_dH'  'MRP41_dV'  'MRP41'
        'MRP42_dH'  'MRP42_dV'  'MRP42'
        'MRP43_dH'  'MRP43_dV'  'MRP43'
        'MRP44_dH'  'MRP44_dV'  'MRP44'
        'MRP45_dH'  'MRP45_dV'  'MRP45'
        'MRP51_dH'  'MRP51_dV'  'MRP51'
        'MRP52_dH'  'MRP52_dV'  'MRP52'
        'MRP53_dH'  'MRP53_dV'  'MRP53'
        'MRP54_dH'  'MRP54_dV'  'MRP54'
        'MRP55_dH'  'MRP55_dV'  'MRP55'
        'MRP56_dH'  'MRP56_dV'  'MRP56'
        'MRP57_dH'  'MRP57_dV'  'MRP57'
        'MRT11_dH'  'MRT11_dV'  'MRT11'
        'MRT12_dH'  'MRT12_dV'  'MRT12'
        'MRT13_dH'  'MRT13_dV'  'MRT13'
        'MRT14_dH'  'MRT14_dV'  'MRT14'
        'MRT15_dH'  'MRT15_dV'  'MRT15'
        'MRT16_dH'  'MRT16_dV'  'MRT16'
        'MRT21_dH'  'MRT21_dV'  'MRT21'
        'MRT22_dH'  'MRT22_dV'  'MRT22'
        'MRT23_dH'  'MRT23_dV'  'MRT23'
        'MRT24_dH'  'MRT24_dV'  'MRT24'
        'MRT25_dH'  'MRT25_dV'  'MRT25'
        'MRT26_dH'  'MRT26_dV'  'MRT26'
        'MRT27_dH'  'MRT27_dV'  'MRT27'
        'MRT31_dH'  'MRT31_dV'  'MRT31'
        'MRT32_dH'  'MRT32_dV'  'MRT32'
        'MRT33_dH'  'MRT33_dV'  'MRT33'
        'MRT34_dH'  'MRT34_dV'  'MRT34'
        'MRT35_dH'  'MRT35_dV'  'MRT35'
        'MRT36_dH'  'MRT36_dV'  'MRT36'
        'MRT37_dH'  'MRT37_dV'  'MRT37'
        'MRT41_dH'  'MRT41_dV'  'MRT41'
        'MRT42_dH'  'MRT42_dV'  'MRT42'
        'MRT43_dH'  'MRT43_dV'  'MRT43'
        'MRT44_dH'  'MRT44_dV'  'MRT44'
        'MRT45_dH'  'MRT45_dV'  'MRT45'
        'MRT46_dH'  'MRT46_dV'  'MRT46'
        'MRT47_dH'  'MRT47_dV'  'MRT47'
        'MRT51_dH'  'MRT51_dV'  'MRT51'
        'MRT52_dH'  'MRT52_dV'  'MRT52'
        'MRT53_dH'  'MRT53_dV'  'MRT53'
        'MRT54_dH'  'MRT54_dV'  'MRT54'
        'MRT55_dH'  'MRT55_dV'  'MRT55'
        'MRT56_dH'  'MRT56_dV'  'MRT56'
        'MRT57_dH'  'MRT57_dV'  'MRT57'
        'MZC01_dH'  'MZC01_dV'  'MZC01'
        'MZC02_dH'  'MZC02_dV'  'MZC02'
        'MZC03_dH'  'MZC03_dV'  'MZC03'
        'MZC04_dH'  'MZC04_dV'  'MZC04'
        'MZF01_dH'  'MZF01_dV'  'MZF01'
        'MZF02_dH'  'MZF02_dV'  'MZF02'
        'MZF03_dH'  'MZF03_dV'  'MZF03'
        'MZO01_dH'  'MZO01_dV'  'MZO01'
        'MZO02_dH'  'MZO02_dV'  'MZO02'
        'MZO03_dH'  'MZO03_dV'  'MZO03'
        'MZP01_dH'  'MZP01_dV'  'MZP01'
        };
      ctf275_planar_combined = label(:,3);
      label = label(:,1:2);

    case {'neuromag122' 'neuromag122alt'}
      % this is the combination of the two versions (with and without space)
      label = {
        'MEG 001'  'MEG 002'  'MEG 001+002'
        'MEG 003'  'MEG 004'  'MEG 003+004'
        'MEG 005'  'MEG 006'  'MEG 005+006'
        'MEG 007'  'MEG 008'  'MEG 007+008'
        'MEG 009'  'MEG 010'  'MEG 009+010'
        'MEG 011'  'MEG 012'  'MEG 011+012'
        'MEG 013'  'MEG 014'  'MEG 013+014'
        'MEG 015'  'MEG 016'  'MEG 015+016'
        'MEG 017'  'MEG 018'  'MEG 017+018'
        'MEG 019'  'MEG 020'  'MEG 019+020'
        'MEG 021'  'MEG 022'  'MEG 021+022'
        'MEG 023'  'MEG 024'  'MEG 023+024'
        'MEG 025'  'MEG 026'  'MEG 025+026'
        'MEG 027'  'MEG 028'  'MEG 027+028'
        'MEG 029'  'MEG 030'  'MEG 029+030'
        'MEG 031'  'MEG 032'  'MEG 031+032'
        'MEG 033'  'MEG 034'  'MEG 033+034'
        'MEG 035'  'MEG 036'  'MEG 035+036'
        'MEG 037'  'MEG 038'  'MEG 037+038'
        'MEG 039'  'MEG 040'  'MEG 039+040'
        'MEG 041'  'MEG 042'  'MEG 041+042'
        'MEG 043'  'MEG 044'  'MEG 043+044'
        'MEG 045'  'MEG 046'  'MEG 045+046'
        'MEG 047'  'MEG 048'  'MEG 047+048'
        'MEG 049'  'MEG 050'  'MEG 049+050'
        'MEG 051'  'MEG 052'  'MEG 051+052'
        'MEG 053'  'MEG 054'  'MEG 053+054'
        'MEG 055'  'MEG 056'  'MEG 055+056'
        'MEG 057'  'MEG 058'  'MEG 057+058'
        'MEG 059'  'MEG 060'  'MEG 059+060'
        'MEG 061'  'MEG 062'  'MEG 061+062'
        'MEG 063'  'MEG 064'  'MEG 063+064'
        'MEG 065'  'MEG 066'  'MEG 065+066'
        'MEG 067'  'MEG 068'  'MEG 067+068'
        'MEG 069'  'MEG 070'  'MEG 069+070'
        'MEG 071'  'MEG 072'  'MEG 071+072'
        'MEG 073'  'MEG 074'  'MEG 073+074'
        'MEG 075'  'MEG 076'  'MEG 075+076'
        'MEG 077'  'MEG 078'  'MEG 077+078'
        'MEG 079'  'MEG 080'  'MEG 079+080'
        'MEG 081'  'MEG 082'  'MEG 081+082'
        'MEG 083'  'MEG 084'  'MEG 083+084'
        'MEG 085'  'MEG 086'  'MEG 085+086'
        'MEG 087'  'MEG 088'  'MEG 087+088'
        'MEG 089'  'MEG 090'  'MEG 089+090'
        'MEG 091'  'MEG 092'  'MEG 091+092'
        'MEG 093'  'MEG 094'  'MEG 093+094'
        'MEG 095'  'MEG 096'  'MEG 095+096'
        'MEG 097'  'MEG 098'  'MEG 097+098'
        'MEG 099'  'MEG 100'  'MEG 099+100'
        'MEG 101'  'MEG 102'  'MEG 101+102'
        'MEG 103'  'MEG 104'  'MEG 103+104'
        'MEG 105'  'MEG 106'  'MEG 105+106'
        'MEG 107'  'MEG 108'  'MEG 107+108'
        'MEG 109'  'MEG 110'  'MEG 109+110'
        'MEG 111'  'MEG 112'  'MEG 111+112'
        'MEG 113'  'MEG 114'  'MEG 113+114'
        'MEG 115'  'MEG 116'  'MEG 115+116'
        'MEG 117'  'MEG 118'  'MEG 117+118'
        'MEG 119'  'MEG 120'  'MEG 119+120'
        'MEG 121'  'MEG 122'  'MEG 121+122'
        % this is an alternative set of labels without a space in them
        'MEG001'  'MEG002'  'MEG001+002'
        'MEG003'  'MEG004'  'MEG003+004'
        'MEG005'  'MEG006'  'MEG005+006'
        'MEG007'  'MEG008'  'MEG007+008'
        'MEG009'  'MEG010'  'MEG009+010'
        'MEG011'  'MEG012'  'MEG011+012'
        'MEG013'  'MEG014'  'MEG013+014'
        'MEG015'  'MEG016'  'MEG015+016'
        'MEG017'  'MEG018'  'MEG017+018'
        'MEG019'  'MEG020'  'MEG019+020'
        'MEG021'  'MEG022'  'MEG021+022'
        'MEG023'  'MEG024'  'MEG023+024'
        'MEG025'  'MEG026'  'MEG025+026'
        'MEG027'  'MEG028'  'MEG027+028'
        'MEG029'  'MEG030'  'MEG029+030'
        'MEG031'  'MEG032'  'MEG031+032'
        'MEG033'  'MEG034'  'MEG033+034'
        'MEG035'  'MEG036'  'MEG035+036'
        'MEG037'  'MEG038'  'MEG037+038'
        'MEG039'  'MEG040'  'MEG039+040'
        'MEG041'  'MEG042'  'MEG041+042'
        'MEG043'  'MEG044'  'MEG043+044'
        'MEG045'  'MEG046'  'MEG045+046'
        'MEG047'  'MEG048'  'MEG047+048'
        'MEG049'  'MEG050'  'MEG049+050'
        'MEG051'  'MEG052'  'MEG051+052'
        'MEG053'  'MEG054'  'MEG053+054'
        'MEG055'  'MEG056'  'MEG055+056'
        'MEG057'  'MEG058'  'MEG057+058'
        'MEG059'  'MEG060'  'MEG059+060'
        'MEG061'  'MEG062'  'MEG061+062'
        'MEG063'  'MEG064'  'MEG063+064'
        'MEG065'  'MEG066'  'MEG065+066'
        'MEG067'  'MEG068'  'MEG067+068'
        'MEG069'  'MEG070'  'MEG069+070'
        'MEG071'  'MEG072'  'MEG071+072'
        'MEG073'  'MEG074'  'MEG073+074'
        'MEG075'  'MEG076'  'MEG075+076'
        'MEG077'  'MEG078'  'MEG077+078'
        'MEG079'  'MEG080'  'MEG079+080'
        'MEG081'  'MEG082'  'MEG081+082'
        'MEG083'  'MEG084'  'MEG083+084'
        'MEG085'  'MEG086'  'MEG085+086'
        'MEG087'  'MEG088'  'MEG087+088'
        'MEG089'  'MEG090'  'MEG089+090'
        'MEG091'  'MEG092'  'MEG091+092'
        'MEG093'  'MEG094'  'MEG093+094'
        'MEG095'  'MEG096'  'MEG095+096'
        'MEG097'  'MEG098'  'MEG097+098'
        'MEG099'  'MEG100'  'MEG099+100'
        'MEG101'  'MEG102'  'MEG101+102'
        'MEG103'  'MEG104'  'MEG103+104'
        'MEG105'  'MEG106'  'MEG105+106'
        'MEG107'  'MEG108'  'MEG107+108'
        'MEG109'  'MEG110'  'MEG109+110'
        'MEG111'  'MEG112'  'MEG111+112'
        'MEG113'  'MEG114'  'MEG113+114'
        'MEG115'  'MEG116'  'MEG115+116'
        'MEG117'  'MEG118'  'MEG117+118'
        'MEG119'  'MEG120'  'MEG119+120'
        'MEG121'  'MEG122'  'MEG121+122'
        };
      neuromag122_combined = label(:,3);
      neuromag122alt_combined = label(:,3);
      label = label(:,1:2);

    case {'neuromag306' 'neuromag306alt'}
      % this is the combination of the two versions (with and without space)
      label = {
        'MEG 0112'  'MEG 0113'  'MEG 0111'  'MEG 0112+0113'
        'MEG 0122'  'MEG 0123'  'MEG 0121'  'MEG 0122+0123'
        'MEG 0132'  'MEG 0133'  'MEG 0131'  'MEG 0132+0133'
        'MEG 0142'  'MEG 0143'  'MEG 0141'  'MEG 0142+0143'
        'MEG 0212'  'MEG 0213'  'MEG 0211'  'MEG 0212+0213'
        'MEG 0222'  'MEG 0223'  'MEG 0221'  'MEG 0222+0223'
        'MEG 0232'  'MEG 0233'  'MEG 0231'  'MEG 0232+0233'
        'MEG 0242'  'MEG 0243'  'MEG 0241'  'MEG 0242+0243'
        'MEG 0312'  'MEG 0313'  'MEG 0311'  'MEG 0312+0313'
        'MEG 0322'  'MEG 0323'  'MEG 0321'  'MEG 0322+0323'
        'MEG 0332'  'MEG 0333'  'MEG 0331'  'MEG 0332+0333'
        'MEG 0342'  'MEG 0343'  'MEG 0341'  'MEG 0342+0343'
        'MEG 0412'  'MEG 0413'  'MEG 0411'  'MEG 0412+0413'
        'MEG 0422'  'MEG 0423'  'MEG 0421'  'MEG 0422+0423'
        'MEG 0432'  'MEG 0433'  'MEG 0431'  'MEG 0432+0433'
        'MEG 0442'  'MEG 0443'  'MEG 0441'  'MEG 0442+0443'
        'MEG 0512'  'MEG 0513'  'MEG 0511'  'MEG 0512+0513'
        'MEG 0522'  'MEG 0523'  'MEG 0521'  'MEG 0522+0523'
        'MEG 0532'  'MEG 0533'  'MEG 0531'  'MEG 0532+0533'
        'MEG 0542'  'MEG 0543'  'MEG 0541'  'MEG 0542+0543'
        'MEG 0612'  'MEG 0613'  'MEG 0611'  'MEG 0612+0613'
        'MEG 0622'  'MEG 0623'  'MEG 0621'  'MEG 0622+0623'
        'MEG 0632'  'MEG 0633'  'MEG 0631'  'MEG 0632+0633'
        'MEG 0642'  'MEG 0643'  'MEG 0641'  'MEG 0642+0643'
        'MEG 0712'  'MEG 0713'  'MEG 0711'  'MEG 0712+0713'
        'MEG 0722'  'MEG 0723'  'MEG 0721'  'MEG 0722+0723'
        'MEG 0732'  'MEG 0733'  'MEG 0731'  'MEG 0732+0733'
        'MEG 0742'  'MEG 0743'  'MEG 0741'  'MEG 0742+0743'
        'MEG 0812'  'MEG 0813'  'MEG 0811'  'MEG 0812+0813'
        'MEG 0822'  'MEG 0823'  'MEG 0821'  'MEG 0822+0823'
        'MEG 0912'  'MEG 0913'  'MEG 0911'  'MEG 0912+0913'
        'MEG 0922'  'MEG 0923'  'MEG 0921'  'MEG 0922+0923'
        'MEG 0932'  'MEG 0933'  'MEG 0931'  'MEG 0932+0933'
        'MEG 0942'  'MEG 0943'  'MEG 0941'  'MEG 0942+0943'
        'MEG 1012'  'MEG 1013'  'MEG 1011'  'MEG 1012+1013'
        'MEG 1022'  'MEG 1023'  'MEG 1021'  'MEG 1022+1023'
        'MEG 1032'  'MEG 1033'  'MEG 1031'  'MEG 1032+1033'
        'MEG 1042'  'MEG 1043'  'MEG 1041'  'MEG 1042+1043'
        'MEG 1112'  'MEG 1113'  'MEG 1111'  'MEG 1112+1113'
        'MEG 1122'  'MEG 1123'  'MEG 1121'  'MEG 1122+1123'
        'MEG 1132'  'MEG 1133'  'MEG 1131'  'MEG 1132+1133'
        'MEG 1142'  'MEG 1143'  'MEG 1141'  'MEG 1142+1143'
        'MEG 1212'  'MEG 1213'  'MEG 1211'  'MEG 1212+1213'
        'MEG 1222'  'MEG 1223'  'MEG 1221'  'MEG 1222+1223'
        'MEG 1232'  'MEG 1233'  'MEG 1231'  'MEG 1232+1233'
        'MEG 1242'  'MEG 1243'  'MEG 1241'  'MEG 1242+1243'
        'MEG 1312'  'MEG 1313'  'MEG 1311'  'MEG 1312+1313'
        'MEG 1322'  'MEG 1323'  'MEG 1321'  'MEG 1322+1323'
        'MEG 1332'  'MEG 1333'  'MEG 1331'  'MEG 1332+1333'
        'MEG 1342'  'MEG 1343'  'MEG 1341'  'MEG 1342+1343'
        'MEG 1412'  'MEG 1413'  'MEG 1411'  'MEG 1412+1413'
        'MEG 1422'  'MEG 1423'  'MEG 1421'  'MEG 1422+1423'
        'MEG 1432'  'MEG 1433'  'MEG 1431'  'MEG 1432+1433'
        'MEG 1442'  'MEG 1443'  'MEG 1441'  'MEG 1442+1443'
        'MEG 1512'  'MEG 1513'  'MEG 1511'  'MEG 1512+1513'
        'MEG 1522'  'MEG 1523'  'MEG 1521'  'MEG 1522+1523'
        'MEG 1532'  'MEG 1533'  'MEG 1531'  'MEG 1532+1533'
        'MEG 1542'  'MEG 1543'  'MEG 1541'  'MEG 1542+1543'
        'MEG 1612'  'MEG 1613'  'MEG 1611'  'MEG 1612+1613'
        'MEG 1622'  'MEG 1623'  'MEG 1621'  'MEG 1622+1623'
        'MEG 1632'  'MEG 1633'  'MEG 1631'  'MEG 1632+1633'
        'MEG 1642'  'MEG 1643'  'MEG 1641'  'MEG 1642+1643'
        'MEG 1712'  'MEG 1713'  'MEG 1711'  'MEG 1712+1713'
        'MEG 1722'  'MEG 1723'  'MEG 1721'  'MEG 1722+1723'
        'MEG 1732'  'MEG 1733'  'MEG 1731'  'MEG 1732+1733'
        'MEG 1742'  'MEG 1743'  'MEG 1741'  'MEG 1742+1743'
        'MEG 1812'  'MEG 1813'  'MEG 1811'  'MEG 1812+1813'
        'MEG 1822'  'MEG 1823'  'MEG 1821'  'MEG 1822+1823'
        'MEG 1832'  'MEG 1833'  'MEG 1831'  'MEG 1832+1833'
        'MEG 1842'  'MEG 1843'  'MEG 1841'  'MEG 1842+1843'
        'MEG 1912'  'MEG 1913'  'MEG 1911'  'MEG 1912+1913'
        'MEG 1922'  'MEG 1923'  'MEG 1921'  'MEG 1922+1923'
        'MEG 1932'  'MEG 1933'  'MEG 1931'  'MEG 1932+1933'
        'MEG 1942'  'MEG 1943'  'MEG 1941'  'MEG 1942+1943'
        'MEG 2012'  'MEG 2013'  'MEG 2011'  'MEG 2012+2013'
        'MEG 2022'  'MEG 2023'  'MEG 2021'  'MEG 2022+2023'
        'MEG 2032'  'MEG 2033'  'MEG 2031'  'MEG 2032+2033'
        'MEG 2042'  'MEG 2043'  'MEG 2041'  'MEG 2042+2043'
        'MEG 2112'  'MEG 2113'  'MEG 2111'  'MEG 2112+2113'
        'MEG 2122'  'MEG 2123'  'MEG 2121'  'MEG 2122+2123'
        'MEG 2132'  'MEG 2133'  'MEG 2131'  'MEG 2132+2133'
        'MEG 2142'  'MEG 2143'  'MEG 2141'  'MEG 2142+2143'
        'MEG 2212'  'MEG 2213'  'MEG 2211'  'MEG 2212+2213'
        'MEG 2222'  'MEG 2223'  'MEG 2221'  'MEG 2222+2223'
        'MEG 2232'  'MEG 2233'  'MEG 2231'  'MEG 2232+2233'
        'MEG 2242'  'MEG 2243'  'MEG 2241'  'MEG 2242+2243'
        'MEG 2312'  'MEG 2313'  'MEG 2311'  'MEG 2312+2313'
        'MEG 2322'  'MEG 2323'  'MEG 2321'  'MEG 2322+2323'
        'MEG 2332'  'MEG 2333'  'MEG 2331'  'MEG 2332+2333'
        'MEG 2342'  'MEG 2343'  'MEG 2341'  'MEG 2342+2343'
        'MEG 2412'  'MEG 2413'  'MEG 2411'  'MEG 2412+2413'
        'MEG 2422'  'MEG 2423'  'MEG 2421'  'MEG 2422+2423'
        'MEG 2432'  'MEG 2433'  'MEG 2431'  'MEG 2432+2433'
        'MEG 2442'  'MEG 2443'  'MEG 2441'  'MEG 2442+2443'
        'MEG 2512'  'MEG 2513'  'MEG 2511'  'MEG 2512+2513'
        'MEG 2522'  'MEG 2523'  'MEG 2521'  'MEG 2522+2523'
        'MEG 2532'  'MEG 2533'  'MEG 2531'  'MEG 2532+2533'
        'MEG 2542'  'MEG 2543'  'MEG 2541'  'MEG 2542+2543'
        'MEG 2612'  'MEG 2613'  'MEG 2611'  'MEG 2612+2613'
        'MEG 2622'  'MEG 2623'  'MEG 2621'  'MEG 2622+2623'
        'MEG 2632'  'MEG 2633'  'MEG 2631'  'MEG 2632+2633'
        'MEG 2642'  'MEG 2643'  'MEG 2641'  'MEG 2642+2643'
        % this is an alternative set of labels without a space in them
        'MEG0112'  'MEG0113'  'MEG0111'  'MEG0112+0113'
        'MEG0122'  'MEG0123'  'MEG0121'  'MEG0122+0123'
        'MEG0132'  'MEG0133'  'MEG0131'  'MEG0132+0133'
        'MEG0142'  'MEG0143'  'MEG0141'  'MEG0142+0143'
        'MEG0212'  'MEG0213'  'MEG0211'  'MEG0212+0213'
        'MEG0222'  'MEG0223'  'MEG0221'  'MEG0222+0223'
        'MEG0232'  'MEG0233'  'MEG0231'  'MEG0232+0233'
        'MEG0242'  'MEG0243'  'MEG0241'  'MEG0242+0243'
        'MEG0312'  'MEG0313'  'MEG0311'  'MEG0312+0313'
        'MEG0322'  'MEG0323'  'MEG0321'  'MEG0322+0323'
        'MEG0332'  'MEG0333'  'MEG0331'  'MEG0332+0333'
        'MEG0342'  'MEG0343'  'MEG0341'  'MEG0342+0343'
        'MEG0412'  'MEG0413'  'MEG0411'  'MEG0412+0413'
        'MEG0422'  'MEG0423'  'MEG0421'  'MEG0422+0423'
        'MEG0432'  'MEG0433'  'MEG0431'  'MEG0432+0433'
        'MEG0442'  'MEG0443'  'MEG0441'  'MEG0442+0443'
        'MEG0512'  'MEG0513'  'MEG0511'  'MEG0512+0513'
        'MEG0522'  'MEG0523'  'MEG0521'  'MEG0522+0523'
        'MEG0532'  'MEG0533'  'MEG0531'  'MEG0532+0533'
        'MEG0542'  'MEG0543'  'MEG0541'  'MEG0542+0543'
        'MEG0612'  'MEG0613'  'MEG0611'  'MEG0612+0613'
        'MEG0622'  'MEG0623'  'MEG0621'  'MEG0622+0623'
        'MEG0632'  'MEG0633'  'MEG0631'  'MEG0632+0633'
        'MEG0642'  'MEG0643'  'MEG0641'  'MEG0642+0643'
        'MEG0712'  'MEG0713'  'MEG0711'  'MEG0712+0713'
        'MEG0722'  'MEG0723'  'MEG0721'  'MEG0722+0723'
        'MEG0732'  'MEG0733'  'MEG0731'  'MEG0732+0733'
        'MEG0742'  'MEG0743'  'MEG0741'  'MEG0742+0743'
        'MEG0812'  'MEG0813'  'MEG0811'  'MEG0812+0813'
        'MEG0822'  'MEG0823'  'MEG0821'  'MEG0822+0823'
        'MEG0912'  'MEG0913'  'MEG0911'  'MEG0912+0913'
        'MEG0922'  'MEG0923'  'MEG0921'  'MEG0922+0923'
        'MEG0932'  'MEG0933'  'MEG0931'  'MEG0932+0933'
        'MEG0942'  'MEG0943'  'MEG0941'  'MEG0942+0943'
        'MEG1012'  'MEG1013'  'MEG1011'  'MEG1012+1013'
        'MEG1022'  'MEG1023'  'MEG1021'  'MEG1022+1023'
        'MEG1032'  'MEG1033'  'MEG1031'  'MEG1032+1033'
        'MEG1042'  'MEG1043'  'MEG1041'  'MEG1042+1043'
        'MEG1112'  'MEG1113'  'MEG1111'  'MEG1112+1113'
        'MEG1122'  'MEG1123'  'MEG1121'  'MEG1122+1123'
        'MEG1132'  'MEG1133'  'MEG1131'  'MEG1132+1133'
        'MEG1142'  'MEG1143'  'MEG1141'  'MEG1142+1143'
        'MEG1212'  'MEG1213'  'MEG1211'  'MEG1212+1213'
        'MEG1222'  'MEG1223'  'MEG1221'  'MEG1222+1223'
        'MEG1232'  'MEG1233'  'MEG1231'  'MEG1232+1233'
        'MEG1242'  'MEG1243'  'MEG1241'  'MEG1242+1243'
        'MEG1312'  'MEG1313'  'MEG1311'  'MEG1312+1313'
        'MEG1322'  'MEG1323'  'MEG1321'  'MEG1322+1323'
        'MEG1332'  'MEG1333'  'MEG1331'  'MEG1332+1333'
        'MEG1342'  'MEG1343'  'MEG1341'  'MEG1342+1343'
        'MEG1412'  'MEG1413'  'MEG1411'  'MEG1412+1413'
        'MEG1422'  'MEG1423'  'MEG1421'  'MEG1422+1423'
        'MEG1432'  'MEG1433'  'MEG1431'  'MEG1432+1433'
        'MEG1442'  'MEG1443'  'MEG1441'  'MEG1442+1443'
        'MEG1512'  'MEG1513'  'MEG1511'  'MEG1512+1513'
        'MEG1522'  'MEG1523'  'MEG1521'  'MEG1522+1523'
        'MEG1532'  'MEG1533'  'MEG1531'  'MEG1532+1533'
        'MEG1542'  'MEG1543'  'MEG1541'  'MEG1542+1543'
        'MEG1612'  'MEG1613'  'MEG1611'  'MEG1612+1613'
        'MEG1622'  'MEG1623'  'MEG1621'  'MEG1622+1623'
        'MEG1632'  'MEG1633'  'MEG1631'  'MEG1632+1633'
        'MEG1642'  'MEG1643'  'MEG1641'  'MEG1642+1643'
        'MEG1712'  'MEG1713'  'MEG1711'  'MEG1712+1713'
        'MEG1722'  'MEG1723'  'MEG1721'  'MEG1722+1723'
        'MEG1732'  'MEG1733'  'MEG1731'  'MEG1732+1733'
        'MEG1742'  'MEG1743'  'MEG1741'  'MEG1742+1743'
        'MEG1812'  'MEG1813'  'MEG1811'  'MEG1812+1813'
        'MEG1822'  'MEG1823'  'MEG1821'  'MEG1822+1823'
        'MEG1832'  'MEG1833'  'MEG1831'  'MEG1832+1833'
        'MEG1842'  'MEG1843'  'MEG1841'  'MEG1842+1843'
        'MEG1912'  'MEG1913'  'MEG1911'  'MEG1912+1913'
        'MEG1922'  'MEG1923'  'MEG1921'  'MEG1922+1923'
        'MEG1932'  'MEG1933'  'MEG1931'  'MEG1932+1933'
        'MEG1942'  'MEG1943'  'MEG1941'  'MEG1942+1943'
        'MEG2012'  'MEG2013'  'MEG2011'  'MEG2012+2013'
        'MEG2022'  'MEG2023'  'MEG2021'  'MEG2022+2023'
        'MEG2032'  'MEG2033'  'MEG2031'  'MEG2032+2033'
        'MEG2042'  'MEG2043'  'MEG2041'  'MEG2042+2043'
        'MEG2112'  'MEG2113'  'MEG2111'  'MEG2112+2113'
        'MEG2122'  'MEG2123'  'MEG2121'  'MEG2122+2123'
        'MEG2132'  'MEG2133'  'MEG2131'  'MEG2132+2133'
        'MEG2142'  'MEG2143'  'MEG2141'  'MEG2142+2143'
        'MEG2212'  'MEG2213'  'MEG2211'  'MEG2212+2213'
        'MEG2222'  'MEG2223'  'MEG2221'  'MEG2222+2223'
        'MEG2232'  'MEG2233'  'MEG2231'  'MEG2232+2233'
        'MEG2242'  'MEG2243'  'MEG2241'  'MEG2242+2243'
        'MEG2312'  'MEG2313'  'MEG2311'  'MEG2312+2313'
        'MEG2322'  'MEG2323'  'MEG2321'  'MEG2322+2323'
        'MEG2332'  'MEG2333'  'MEG2331'  'MEG2332+2333'
        'MEG2342'  'MEG2343'  'MEG2341'  'MEG2342+2343'
        'MEG2412'  'MEG2413'  'MEG2411'  'MEG2412+2413'
        'MEG2422'  'MEG2423'  'MEG2421'  'MEG2422+2423'
        'MEG2432'  'MEG2433'  'MEG2431'  'MEG2432+2433'
        'MEG2442'  'MEG2443'  'MEG2441'  'MEG2442+2443'
        'MEG2512'  'MEG2513'  'MEG2511'  'MEG2512+2513'
        'MEG2522'  'MEG2523'  'MEG2521'  'MEG2522+2523'
        'MEG2532'  'MEG2533'  'MEG2531'  'MEG2532+2533'
        'MEG2542'  'MEG2543'  'MEG2541'  'MEG2542+2543'
        'MEG2612'  'MEG2613'  'MEG2611'  'MEG2612+2613'
        'MEG2622'  'MEG2623'  'MEG2621'  'MEG2622+2623'
        'MEG2632'  'MEG2633'  'MEG2631'  'MEG2632+2633'
        'MEG2642'  'MEG2643'  'MEG2641'  'MEG2642+2643'
        };
      neuromag306_combined = label(:,4);
      neuromag306alt_combined = label(:,4);
      label = label(:,1:3);

    case 'eeg1020'
      label = {
        'Fp1'
        'Fpz'
        'Fp2'
        'F7'
        'F3'
        'Fz'
        'F4'
        'F8'
        'T7'
        'C3'
        'Cz'
        'C4'
        'T8'
        'P7'
        'P3'
        'Pz'
        'P4'
        'P8'
        'O1'
        'Oz'
        'O2'};

      % Add also reference and some alternative labels that might be used
      label = cat(1, label, {'A1' 'A2' 'M1' 'M2' 'T3' 'T4' 'T5' 'T6'}');

    case 'eeg1010'
      label = {
        'Fp1'
        'Fpz'
        'Fp2'
        'AF9'
        'AF7'
        'AF5'
        'AF3'
        'AF1'
        'AFz'
        'AF2'
        'AF4'
        'AF6'
        'AF8'
        'AF10'
        'F9'
        'F7'
        'F5'
        'F3'
        'F1'
        'Fz'
        'F2'
        'F4'
        'F6'
        'F8'
        'F10'
        'FT9'
        'FT7'
        'FC5'
        'FC3'
        'FC1'
        'FCz'
        'FC2'
        'FC4'
        'FC6'
        'FT8'
        'FT10'
        'T9'
        'T7'
        'C5'
        'C3'
        'C1'
        'Cz'
        'C2'
        'C4'
        'C6'
        'T8'
        'T10'
        'TP9'
        'TP7'
        'CP5'
        'CP3'
        'CP1'
        'CPz'
        'CP2'
        'CP4'
        'CP6'
        'TP8'
        'TP10'
        'P9'
        'P7'
        'P5'
        'P3'
        'P1'
        'Pz'
        'P2'
        'P4'
        'P6'
        'P8'
        'P10'
        'PO9'
        'PO7'
        'PO5'
        'PO3'
        'PO1'
        'POz'
        'PO2'
        'PO4'
        'PO6'
        'PO8'
        'PO10'
        'O1'
        'Oz'
        'O2'
        'I1'
        'Iz'
        'I2'
        };

      % Add also reference and some alternative labels that might be used
      label = cat(1, label, {'A1' 'A2' 'M1' 'M2' 'T3' 'T4' 'T5' 'T6'}');

    case 'eeg1005'
      label = {
        'Fp1'
        'Fpz'
        'Fp2'
        'AF9'
        'AF7'
        'AF5'
        'AF3'
        'AF1'
        'AFz'
        'AF2'
        'AF4'
        'AF6'
        'AF8'
        'AF10'
        'F9'
        'F7'
        'F5'
        'F3'
        'F1'
        'Fz'
        'F2'
        'F4'
        'F6'
        'F8'
        'F10'
        'FT9'
        'FT7'
        'FC5'
        'FC3'
        'FC1'
        'FCz'
        'FC2'
        'FC4'
        'FC6'
        'FT8'
        'FT10'
        'T9'
        'T7'
        'C5'
        'C3'
        'C1'
        'Cz'
        'C2'
        'C4'
        'C6'
        'T8'
        'T10'
        'TP9'
        'TP7'
        'CP5'
        'CP3'
        'CP1'
        'CPz'
        'CP2'
        'CP4'
        'CP6'
        'TP8'
        'TP10'
        'P9'
        'P7'
        'P5'
        'P3'
        'P1'
        'Pz'
        'P2'
        'P4'
        'P6'
        'P8'
        'P10'
        'PO9'
        'PO7'
        'PO5'
        'PO3'
        'PO1'
        'POz'
        'PO2'
        'PO4'
        'PO6'
        'PO8'
        'PO10'
        'O1'
        'Oz'
        'O2'
        'I1'
        'Iz'
        'I2'
        'AFp9h'
        'AFp7h'
        'AFp5h'
        'AFp3h'
        'AFp1h'
        'AFp2h'
        'AFp4h'
        'AFp6h'
        'AFp8h'
        'AFp10h'
        'AFF9h'
        'AFF7h'
        'AFF5h'
        'AFF3h'
        'AFF1h'
        'AFF2h'
        'AFF4h'
        'AFF6h'
        'AFF8h'
        'AFF10h'
        'FFT9h'
        'FFT7h'
        'FFC5h'
        'FFC3h'
        'FFC1h'
        'FFC2h'
        'FFC4h'
        'FFC6h'
        'FFT8h'
        'FFT10h'
        'FTT9h'
        'FTT7h'
        'FCC5h'
        'FCC3h'
        'FCC1h'
        'FCC2h'
        'FCC4h'
        'FCC6h'
        'FTT8h'
        'FTT10h'
        'TTP9h'
        'TTP7h'
        'CCP5h'
        'CCP3h'
        'CCP1h'
        'CCP2h'
        'CCP4h'
        'CCP6h'
        'TTP8h'
        'TTP10h'
        'TPP9h'
        'TPP7h'
        'CPP5h'
        'CPP3h'
        'CPP1h'
        'CPP2h'
        'CPP4h'
        'CPP6h'
        'TPP8h'
        'TPP10h'
        'PPO9h'
        'PPO7h'
        'PPO5h'
        'PPO3h'
        'PPO1h'
        'PPO2h'
        'PPO4h'
        'PPO6h'
        'PPO8h'
        'PPO10h'
        'POO9h'
        'POO7h'
        'POO5h'
        'POO3h'
        'POO1h'
        'POO2h'
        'POO4h'
        'POO6h'
        'POO8h'
        'POO10h'
        'OI1h'
        'OI2h'
        'Fp1h'
        'Fp2h'
        'AF9h'
        'AF7h'
        'AF5h'
        'AF3h'
        'AF1h'
        'AF2h'
        'AF4h'
        'AF6h'
        'AF8h'
        'AF10h'
        'F9h'
        'F7h'
        'F5h'
        'F3h'
        'F1h'
        'F2h'
        'F4h'
        'F6h'
        'F8h'
        'F10h'
        'FT9h'
        'FT7h'
        'FC5h'
        'FC3h'
        'FC1h'
        'FC2h'
        'FC4h'
        'FC6h'
        'FT8h'
        'FT10h'
        'T9h'
        'T7h'
        'C5h'
        'C3h'
        'C1h'
        'C2h'
        'C4h'
        'C6h'
        'T8h'
        'T10h'
        'TP9h'
        'TP7h'
        'CP5h'
        'CP3h'
        'CP1h'
        'CP2h'
        'CP4h'
        'CP6h'
        'TP8h'
        'TP10h'
        'P9h'
        'P7h'
        'P5h'
        'P3h'
        'P1h'
        'P2h'
        'P4h'
        'P6h'
        'P8h'
        'P10h'
        'PO9h'
        'PO7h'
        'PO5h'
        'PO3h'
        'PO1h'
        'PO2h'
        'PO4h'
        'PO6h'
        'PO8h'
        'PO10h'
        'O1h'
        'O2h'
        'I1h'
        'I2h'
        'AFp9'
        'AFp7'
        'AFp5'
        'AFp3'
        'AFp1'
        'AFpz'
        'AFp2'
        'AFp4'
        'AFp6'
        'AFp8'
        'AFp10'
        'AFF9'
        'AFF7'
        'AFF5'
        'AFF3'
        'AFF1'
        'AFFz'
        'AFF2'
        'AFF4'
        'AFF6'
        'AFF8'
        'AFF10'
        'FFT9'
        'FFT7'
        'FFC5'
        'FFC3'
        'FFC1'
        'FFCz'
        'FFC2'
        'FFC4'
        'FFC6'
        'FFT8'
        'FFT10'
        'FTT9'
        'FTT7'
        'FCC5'
        'FCC3'
        'FCC1'
        'FCCz'
        'FCC2'
        'FCC4'
        'FCC6'
        'FTT8'
        'FTT10'
        'TTP9'
        'TTP7'
        'CCP5'
        'CCP3'
        'CCP1'
        'CCPz'
        'CCP2'
        'CCP4'
        'CCP6'
        'TTP8'
        'TTP10'
        'TPP9'
        'TPP7'
        'CPP5'
        'CPP3'
        'CPP1'
        'CPPz'
        'CPP2'
        'CPP4'
        'CPP6'
        'TPP8'
        'TPP10'
        'PPO9'
        'PPO7'
        'PPO5'
        'PPO3'
        'PPO1'
        'PPOz'
        'PPO2'
        'PPO4'
        'PPO6'
        'PPO8'
        'PPO10'
        'POO9'
        'POO7'
        'POO5'
        'POO3'
        'POO1'
        'POOz'
        'POO2'
        'POO4'
        'POO6'
        'POO8'
        'POO10'
        'OI1'
        'OIz'
        'OI2'
        };

      % Add also reference and some alternative labels that might be used
      label = cat(1, label, {'A1' 'A2' 'M1' 'M2' 'T3' 'T4' 'T5' 'T6'}');

    case 'ext1020'
      % start with the eeg1005 list
      label = {
        'Fp1'
        'Fpz'
        'Fp2'
        'AF9'
        'AF7'
        'AF5'
        'AF3'
        'AF1'
        'AFz'
        'AF2'
        'AF4'
        'AF6'
        'AF8'
        'AF10'
        'F9'
        'F7'
        'F5'
        'F3'
        'F1'
        'Fz'
        'F2'
        'F4'
        'F6'
        'F8'
        'F10'
        'FT9'
        'FT7'
        'FC5'
        'FC3'
        'FC1'
        'FCz'
        'FC2'
        'FC4'
        'FC6'
        'FT8'
        'FT10'
        'T9'
        'T7'
        'C5'
        'C3'
        'C1'
        'Cz'
        'C2'
        'C4'
        'C6'
        'T8'
        'T10'
        'TP9'
        'TP7'
        'CP5'
        'CP3'
        'CP1'
        'CPz'
        'CP2'
        'CP4'
        'CP6'
        'TP8'
        'TP10'
        'P9'
        'P7'
        'P5'
        'P3'
        'P1'
        'Pz'
        'P2'
        'P4'
        'P6'
        'P8'
        'P10'
        'PO9'
        'PO7'
        'PO5'
        'PO3'
        'PO1'
        'POz'
        'PO2'
        'PO4'
        'PO6'
        'PO8'
        'PO10'
        'O1'
        'Oz'
        'O2'
        'I1'
        'Iz'
        'I2'
        'AFp9h'
        'AFp7h'
        'AFp5h'
        'AFp3h'
        'AFp1h'
        'AFp2h'
        'AFp4h'
        'AFp6h'
        'AFp8h'
        'AFp10h'
        'AFF9h'
        'AFF7h'
        'AFF5h'
        'AFF3h'
        'AFF1h'
        'AFF2h'
        'AFF4h'
        'AFF6h'
        'AFF8h'
        'AFF10h'
        'FFT9h'
        'FFT7h'
        'FFC5h'
        'FFC3h'
        'FFC1h'
        'FFC2h'
        'FFC4h'
        'FFC6h'
        'FFT8h'
        'FFT10h'
        'FTT9h'
        'FTT7h'
        'FCC5h'
        'FCC3h'
        'FCC1h'
        'FCC2h'
        'FCC4h'
        'FCC6h'
        'FTT8h'
        'FTT10h'
        'TTP9h'
        'TTP7h'
        'CCP5h'
        'CCP3h'
        'CCP1h'
        'CCP2h'
        'CCP4h'
        'CCP6h'
        'TTP8h'
        'TTP10h'
        'TPP9h'
        'TPP7h'
        'CPP5h'
        'CPP3h'
        'CPP1h'
        'CPP2h'
        'CPP4h'
        'CPP6h'
        'TPP8h'
        'TPP10h'
        'PPO9h'
        'PPO7h'
        'PPO5h'
        'PPO3h'
        'PPO1h'
        'PPO2h'
        'PPO4h'
        'PPO6h'
        'PPO8h'
        'PPO10h'
        'POO9h'
        'POO7h'
        'POO5h'
        'POO3h'
        'POO1h'
        'POO2h'
        'POO4h'
        'POO6h'
        'POO8h'
        'POO10h'
        'OI1h'
        'OI2h'
        'Fp1h'
        'Fp2h'
        'AF9h'
        'AF7h'
        'AF5h'
        'AF3h'
        'AF1h'
        'AF2h'
        'AF4h'
        'AF6h'
        'AF8h'
        'AF10h'
        'F9h'
        'F7h'
        'F5h'
        'F3h'
        'F1h'
        'F2h'
        'F4h'
        'F6h'
        'F8h'
        'F10h'
        'FT9h'
        'FT7h'
        'FC5h'
        'FC3h'
        'FC1h'
        'FC2h'
        'FC4h'
        'FC6h'
        'FT8h'
        'FT10h'
        'T9h'
        'T7h'
        'C5h'
        'C3h'
        'C1h'
        'C2h'
        'C4h'
        'C6h'
        'T8h'
        'T10h'
        'TP9h'
        'TP7h'
        'CP5h'
        'CP3h'
        'CP1h'
        'CP2h'
        'CP4h'
        'CP6h'
        'TP8h'
        'TP10h'
        'P9h'
        'P7h'
        'P5h'
        'P3h'
        'P1h'
        'P2h'
        'P4h'
        'P6h'
        'P8h'
        'P10h'
        'PO9h'
        'PO7h'
        'PO5h'
        'PO3h'
        'PO1h'
        'PO2h'
        'PO4h'
        'PO6h'
        'PO8h'
        'PO10h'
        'O1h'
        'O2h'
        'I1h'
        'I2h'
        'AFp9'
        'AFp7'
        'AFp5'
        'AFp3'
        'AFp1'
        'AFpz'
        'AFp2'
        'AFp4'
        'AFp6'
        'AFp8'
        'AFp10'
        'AFF9'
        'AFF7'
        'AFF5'
        'AFF3'
        'AFF1'
        'AFFz'
        'AFF2'
        'AFF4'
        'AFF6'
        'AFF8'
        'AFF10'
        'FFT9'
        'FFT7'
        'FFC5'
        'FFC3'
        'FFC1'
        'FFCz'
        'FFC2'
        'FFC4'
        'FFC6'
        'FFT8'
        'FFT10'
        'FTT9'
        'FTT7'
        'FCC5'
        'FCC3'
        'FCC1'
        'FCCz'
        'FCC2'
        'FCC4'
        'FCC6'
        'FTT8'
        'FTT10'
        'TTP9'
        'TTP7'
        'CCP5'
        'CCP3'
        'CCP1'
        'CCPz'
        'CCP2'
        'CCP4'
        'CCP6'
        'TTP8'
        'TTP10'
        'TPP9'
        'TPP7'
        'CPP5'
        'CPP3'
        'CPP1'
        'CPPz'
        'CPP2'
        'CPP4'
        'CPP6'
        'TPP8'
        'TPP10'
        'PPO9'
        'PPO7'
        'PPO5'
        'PPO3'
        'PPO1'
        'PPOz'
        'PPO2'
        'PPO4'
        'PPO6'
        'PPO8'
        'PPO10'
        'POO9'
        'POO7'
        'POO5'
        'POO3'
        'POO1'
        'POOz'
        'POO2'
        'POO4'
        'POO6'
        'POO8'
        'POO10'
        'OI1'
        'OIz'
        'OI2'
        };

      % Add also reference and some alternative labels that might be used
      label = cat(1, label, {'A1' 'A2' 'M1' 'M2' 'T3' 'T4' 'T5' 'T6'}');

      % This is to account for all variants of case in 1020 systems
      label = unique(cat(1, label, upper(label), lower(label)));

    case 'biosemi64'
      label = {
        'A1'
        'A2'
        'A3'
        'A4'
        'A5'
        'A6'
        'A7'
        'A8'
        'A9'
        'A10'
        'A11'
        'A12'
        'A13'
        'A14'
        'A15'
        'A16'
        'A17'
        'A18'
        'A19'
        'A20'
        'A21'
        'A22'
        'A23'
        'A24'
        'A25'
        'A26'
        'A27'
        'A28'
        'A29'
        'A30'
        'A31'
        'A32'
        'B1'
        'B2'
        'B3'
        'B4'
        'B5'
        'B6'
        'B7'
        'B8'
        'B9'
        'B10'
        'B11'
        'B12'
        'B13'
        'B14'
        'B15'
        'B16'
        'B17'
        'B18'
        'B19'
        'B20'
        'B21'
        'B22'
        'B23'
        'B24'
        'B25'
        'B26'
        'B27'
        'B28'
        'B29'
        'B30'
        'B31'
        'B32'
        };

    case 'biosemi128'
      label = {
        'A1'
        'A2'
        'A3'
        'A4'
        'A5'
        'A6'
        'A7'
        'A8'
        'A9'
        'A10'
        'A11'
        'A12'
        'A13'
        'A14'
        'A15'
        'A16'
        'A17'
        'A18'
        'A19'
        'A20'
        'A21'
        'A22'
        'A23'
        'A24'
        'A25'
        'A26'
        'A27'
        'A28'
        'A29'
        'A30'
        'A31'
        'A32'
        'B1'
        'B2'
        'B3'
        'B4'
        'B5'
        'B6'
        'B7'
        'B8'
        'B9'
        'B10'
        'B11'
        'B12'
        'B13'
        'B14'
        'B15'
        'B16'
        'B17'
        'B18'
        'B19'
        'B20'
        'B21'
        'B22'
        'B23'
        'B24'
        'B25'
        'B26'
        'B27'
        'B28'
        'B29'
        'B30'
        'B31'
        'B32'
        'C1'
        'C2'
        'C3'
        'C4'
        'C5'
        'C6'
        'C7'
        'C8'
        'C9'
        'C10'
        'C11'
        'C12'
        'C13'
        'C14'
        'C15'
        'C16'
        'C17'
        'C18'
        'C19'
        'C20'
        'C21'
        'C22'
        'C23'
        'C24'
        'C25'
        'C26'
        'C27'
        'C28'
        'C29'
        'C30'
        'C31'
        'C32'
        'D1'
        'D2'
        'D3'
        'D4'
        'D5'
        'D6'
        'D7'
        'D8'
        'D9'
        'D10'
        'D11'
        'D12'
        'D13'
        'D14'
        'D15'
        'D16'
        'D17'
        'D18'
        'D19'
        'D20'
        'D21'
        'D22'
        'D23'
        'D24'
        'D25'
        'D26'
        'D27'
        'D28'
        'D29'
        'D30'
        'D31'
        'D32'
        };

    case 'biosemi256'
      label = {
        'A1'
        'A2'
        'A3'
        'A4'
        'A5'
        'A6'
        'A7'
        'A8'
        'A9'
        'A10'
        'A11'
        'A12'
        'A13'
        'A14'
        'A15'
        'A16'
        'A17'
        'A18'
        'A19'
        'A20'
        'A21'
        'A22'
        'A23'
        'A24'
        'A25'
        'A26'
        'A27'
        'A28'
        'A29'
        'A30'
        'A31'
        'A32'
        'B1'
        'B2'
        'B3'
        'B4'
        'B5'
        'B6'
        'B7'
        'B8'
        'B9'
        'B10'
        'B11'
        'B12'
        'B13'
        'B14'
        'B15'
        'B16'
        'B17'
        'B18'
        'B19'
        'B20'
        'B21'
        'B22'
        'B23'
        'B24'
        'B25'
        'B26'
        'B27'
        'B28'
        'B29'
        'B30'
        'B31'
        'B32'
        'C1'
        'C2'
        'C3'
        'C4'
        'C5'
        'C6'
        'C7'
        'C8'
        'C9'
        'C10'
        'C11'
        'C12'
        'C13'
        'C14'
        'C15'
        'C16'
        'C17'
        'C18'
        'C19'
        'C20'
        'C21'
        'C22'
        'C23'
        'C24'
        'C25'
        'C26'
        'C27'
        'C28'
        'C29'
        'C30'
        'C31'
        'C32'
        'D1'
        'D2'
        'D3'
        'D4'
        'D5'
        'D6'
        'D7'
        'D8'
        'D9'
        'D10'
        'D11'
        'D12'
        'D13'
        'D14'
        'D15'
        'D16'
        'D17'
        'D18'
        'D19'
        'D20'
        'D21'
        'D22'
        'D23'
        'D24'
        'D25'
        'D26'
        'D27'
        'D28'
        'D29'
        'D30'
        'D31'
        'D32'
        'E1'
        'E2'
        'E3'
        'E4'
        'E5'
        'E6'
        'E7'
        'E8'
        'E9'
        'E10'
        'E11'
        'E12'
        'E13'
        'E14'
        'E15'
        'E16'
        'E17'
        'E18'
        'E19'
        'E20'
        'E21'
        'E22'
        'E23'
        'E24'
        'E25'
        'E26'
        'E27'
        'E28'
        'E29'
        'E30'
        'E31'
        'E32'
        'F1'
        'F2'
        'F3'
        'F4'
        'F5'
        'F6'
        'F7'
        'F8'
        'F9'
        'F10'
        'F11'
        'F12'
        'F13'
        'F14'
        'F15'
        'F16'
        'F17'
        'F18'
        'F19'
        'F20'
        'F21'
        'F22'
        'F23'
        'F24'
        'F25'
        'F26'
        'F27'
        'F28'
        'F29'
        'F30'
        'F31'
        'F32'
        'G1'
        'G2'
        'G3'
        'G4'
        'G5'
        'G6'
        'G7'
        'G8'
        'G9'
        'G10'
        'G11'
        'G12'
        'G13'
        'G14'
        'G15'
        'G16'
        'G17'
        'G18'
        'G19'
        'G20'
        'G21'
        'G22'
        'G23'
        'G24'
        'G25'
        'G26'
        'G27'
        'G28'
        'G29'
        'G30'
        'G31'
        'G32'
        'H1'
        'H2'
        'H3'
        'H4'
        'H5'
        'H6'
        'H7'
        'H8'
        'H9'
        'H10'
        'H11'
        'H12'
        'H13'
        'H14'
        'H15'
        'H16'
        'H17'
        'H18'
        'H19'
        'H20'
        'H21'
        'H22'
        'H23'
        'H24'
        'H25'
        'H26'
        'H27'
        'H28'
        'H29'
        'H30'
        'H31'
        'H32'
        };

    case 'egi32'
      % this should be  uppercase for consistency with ft_read_header
      label = cell(33, 1);
      for i = 1:33
        label{i} = sprintf('E%d', i);
      end
      % there might also be a reference channel, but its name is inconsistent
      % it might be Cz, REF, VREF or 'vertex reference'

    case 'egi64'
      % this should be  uppercase for consistency with ft_read_header
      label = cell(65, 1);
      for i = 1:65
        label{i} = sprintf('E%d', i);
      end
      % there might also be a reference channel, but its name is inconsistent
      % it might be Cz, REF, VREF or 'vertex reference'

    case 'egi128'
      % this should be  uppercase for consistency with ft_read_header
      label = cell(129, 1);
      for i = 1:129
        label{i} = sprintf('E%d', i);
      end
      % there might also be a reference channel, but its name is inconsistent
      % it might be Cz, REF, VREF or 'vertex reference'

    case 'egi256'
      % this should be  uppercase for consistency with ft_read_header
      label = cell(257, 1);
      for i = 1:257
        label{i} = sprintf('E%d', i);
      end
      % there might also be a reference channel, but its name is inconsistent
      % it might be Cz, REF, VREF or 'vertex reference'

    case 'itab28'
      label = {
        'MAG_1'
        'MAG_2'
        'MAG_3'
        'MAG_5'
        'MAG_7'
        'MAG_8'
        'MAG_9'
        'MAG_11'
        'MAG_12'
        'MAG_13'
        'MAG_15'
        'MAG_17'
        'MAG_18'
        'MAG_21'
        'MAG_22'
        'MAG_23'
        'MAG_25'
        'MAG_26'
        'MAG_27'
        'MAG_28'
        };

    case 'itab153'
      label = cell(153,1);
      for i=1:153
        % channel names start counting at zero
        label{i} = sprintf('MAG_%03d',  i-1);
      end

    case 'itab153_planar'
      label = cell(153,3);
      for i=1:153
        % channel names start counting at zero
        label{i,1} = sprintf('MAG_%03d_dH', i-1);
        label{i,2} = sprintf('MAG_%03d_dV', i-1);
        label{i,3} = sprintf('MAG_%03d',  i-1);
      end
      itab153_planar_combined = label(:,3);
      label = label(:,1:2);

    case 'yokogawa9'
      % note that this uses MATLAB style 1-offset indexing and not C style 0-offset indexing
      % this should be consistent with: read_yokogawa_header, ft_channelselection, yokogawa2grad
      label = cell(9,1);
      for i=1:9
        label{i} = sprintf('M%03d',  i);
      end

    case 'yokogawa64'
      % note that this uses MATLAB style 1-offset indexing and not C style 0-offset indexing
      % this should be consistent with: read_yokogawa_header, ft_channelselection, yokogawa2grad
      label = cell(64,1);
      for i=1:64
        label{i} = sprintf('AG%03d', i);
      end

    case 'yokogawa64_planar'
      % note that this uses MATLAB style 1-offset indexing and not C style 0-offset indexing
      % this should be consistent with: read_yokogawa_header, ft_channelselection, yokogawa2grad
      label = cell(64,3);
      for i=1:64
        label{i,1} = sprintf('AG%03d_dH', i);
        label{i,2} = sprintf('AG%03d_dV', i);
        label{i,3} = sprintf('AG%03d', i);
      end
      yokogawa64_planar_combined = label(:,3);
      label = label(:,1:2);

    case 'yokogawa160'
      % note that this uses MATLAB style 1-offset indexing and not C style 0-offset indexing
      % this should be consistent with: read_yokogawa_header, ft_channelselection, yokogawa2grad
      label = cell(160,1);
      for i=1:160
        label{i} = sprintf('AG%03d', i);
      end

    case 'yokogawa160_planar'
      % note that this uses MATLAB style 1-offset indexing and not C style 0-offset indexing
      % this should be consistent with: read_yokogawa_header, ft_channelselection, yokogawa2grad
      label = cell(160,2);
      for i=1:160
        label{i,1} = sprintf('AG%03d_dH', i);
        label{i,2} = sprintf('AG%03d_dV', i);
        label{i,3} = sprintf('AG%03d', i);
      end
      yokogawa160_planar_combined = label(:,3);
      label = label(:,1:2);

    case 'yokogawa440'
      % this should be consistent with read_yokogawa_header, with ft_channelselection and with yokogawa2grad
      label = {
        'AG001'
        'AG002'
        'AG003'
        'AG004'
        'AG005'
        'AG006'
        'AG007'
        'AG008'
        'AG009'
        'AG010'
        'AG011'
        'AG012'
        'AG013'
        'AG014'
        'AG015'
        'AG016'
        'AG017'
        'AG018'
        'AG019'
        'AG020'
        'AG021'
        'AG022'
        'AG023'
        'AG024'
        'AG025'
        'AG026'
        'AG027'
        'AG028'
        'AG029'
        'AG030'
        'AG031'
        'AG032'
        'PG033'
        'PG034'
        'PG035'
        'PG036'
        'PG037'
        'PG038'
        'PG039'
        'PG040'
        'PG041'
        'PG042'
        'PG043'
        'PG044'
        'PG045'
        'PG046'
        'PG047'
        'PG048'
        'PG049'
        'PG050'
        'PG051'
        'PG052'
        'PG053'
        'PG054'
        'PG055'
        'PG056'
        'PG057'
        'PG058'
        'PG059'
        'PG060'
        'PG061'
        'PG062'
        'PG063'
        'PG064'
        'AG065'
        'AG066'
        'AG067'
        'AG068'
        'AG069'
        'AG070'
        'AG071'
        'AG072'
        'AG073'
        'AG074'
        'AG075'
        'AG076'
        'AG077'
        'AG078'
        'AG079'
        'AG080'
        'AG081'
        'AG082'
        'AG083'
        'AG084'
        'AG085'
        'AG086'
        'AG087'
        'AG088'
        'AG089'
        'AG090'
        'AG091'
        'AG092'
        'AG093'
        'AG094'
        'AG095'
        'AG096'
        'PG097'
        'PG098'
        'PG099'
        'PG100'
        'PG101'
        'PG102'
        'PG103'
        'PG104'
        'PG105'
        'PG106'
        'PG107'
        'PG108'
        'PG109'
        'PG110'
        'PG111'
        'PG112'
        'PG113'
        'PG114'
        'PG115'
        'PG116'
        'PG117'
        'PG118'
        'PG119'
        'PG120'
        'PG121'
        'AG122'
        'PG123'
        'PG124'
        'PG125'
        'PG126'
        'PG127'
        'PG128'
        'AG129'
        'AG130'
        'AG131'
        'AG132'
        'AG133'
        'AG134'
        'AG135'
        'AG136'
        'AG137'
        'AG138'
        'AG139'
        'AG140'
        'AG141'
        'AG142'
        'AG143'
        'AG144'
        'AG145'
        'AG146'
        'AG147'
        'AG148'
        'AG149'
        'AG150'
        'AG151'
        'AG152'
        'AG153'
        'AG154'
        'AG155'
        'AG156'
        'AG157'
        'AG158'
        'AG159'
        'AG160'
        'AG161'
        'AG162'
        'AG163'
        'AG164'
        'AG165'
        'PG166'
        'PG167'
        'PG168'
        'PG169'
        'PG170'
        'PG171'
        'PG172'
        'PG173'
        'PG174'
        'PG175'
        'PG176'
        'PG177'
        'PG178'
        'PG179'
        'PG180'
        'PG181'
        'PG182'
        'PG183'
        'PG184'
        'PG185'
        'PG186'
        'PG187'
        'PG188'
        'PG189'
        'PG190'
        'PG191'
        'PG192'
        'AG193'
        'AG194'
        'AG195'
        'AG196'
        'AG197'
        'AG198'
        'AG199'
        'AG200'
        'AG201'
        'AG202'
        'AG203'
        'AG204'
        'AG205'
        'AG206'
        'AG207'
        'AG208'
        'AG209'
        'AG210'
        'AG211'
        'AG212'
        'AG213'
        'AG214'
        'AG215'
        'AG216'
        'AG217'
        'AG218'
        'AG219'
        'AG220'
        'AG221'
        'AG222'
        'AG223'
        'AG224'
        'AG225'
        'AG226'
        'AG227'
        'PG228'
        'PG229'
        'PG230'
        'PG231'
        'PG232'
        'PG233'
        'PG234'
        'PG235'
        'PG236'
        'PG237'
        'PG238'
        'PG239'
        'PG240'
        'PG241'
        'PG242'
        'PG243'
        'PG244'
        'PG245'
        'PG246'
        'PG247'
        'PG248'
        'PG249'
        'PG250'
        'PG251'
        'PG252'
        'PG253'
        'PG254'
        'PG255'
        'PG256'
        'AG257'
        'AG258'
        'AG259'
        'AG260'
        'AG261'
        'AG262'
        'AG263'
        'AG264'
        'AG265'
        'AG266'
        'AG267'
        'AG268'
        'AG269'
        'AG270'
        'AG271'
        'AG272'
        'AG273'
        'AG274'
        'AG275'
        'AG276'
        'AG277'
        'AG278'
        'AG279'
        'AG280'
        'AG281'
        'AG282'
        'AG283'
        'AG284'
        'AG285'
        'AG286'
        'AG287'
        'AG288'
        'PG289'
        'PG290'
        'PG291'
        'PG292'
        'PG293'
        'PG294'
        'PG295'
        'PG296'
        'PG297'
        'PG298'
        'PG299'
        'PG300'
        'PG301'
        'PG302'
        'PG303'
        'PG304'
        'PG305'
        'PG306'
        'PG307'
        'PG308'
        'PG309'
        'PG310'
        'PG311'
        'PG312'
        'PG313'
        'PG314'
        'PG315'
        'PG316'
        'PG317'
        'PG318'
        'PG319'
        'PG320'
        'AG321'
        'AG322'
        'AG323'
        'AG324'
        'AG325'
        'AG326'
        'AG327'
        'AG328'
        'AG329'
        'AG330'
        'AG331'
        'AG332'
        'AG333'
        'AG334'
        'AG335'
        'AG336'
        'AG337'
        'AG338'
        'AG339'
        'AG340'
        'AG341'
        'AG342'
        'AG343'
        'AG344'
        'AG345'
        'AG346'
        'AG347'
        'AG348'
        'AG349'
        'AG350'
        'AG351'
        'AG352'
        'PG353'
        'PG354'
        'PG355'
        'PG356'
        'PG357'
        'PG358'
        'PG359'
        'PG360'
        'PG361'
        'PG362'
        'PG363'
        'PG364'
        'PG365'
        'PG366'
        'PG367'
        'PG368'
        'PG369'
        'PG370'
        'PG371'
        'PG372'
        'PG373'
        'PG374'
        'PG375'
        'PG376'
        'PG377'
        'AG378'
        'PG379'
        'PG380'
        'PG381'
        'PG382'
        'PG383'
        'PG384'
        'AG385'
        'AG386'
        'AG387'
        'AG388'
        'AG389'
        'AG390'
        'AG391'
        'AG392'
        'PG393'
        'PG394'
        'PG395'
        'PG396'
        'PG397'
        'PG398'
        'PG399'
        'PG400'
        'RM401'
        'RM402'
        'RM403'
        'RM404'
        'RM405'
        'RM406'
        'RM407'
        'RM408'
        'RM409'
        'RM410'
        'RM411'
        'RM412'
        };

    case 'yokogawa440_planar'
      % this should be consistent with read_yokogawa_header, with
      % ft_channelselection and with yokogawa2grad
      label = {
        'AG001_dH'  'AG001_dV'  'AG001'
        'AG002_dH'  'AG002_dV'  'AG002'
        'AG003_dH'  'AG003_dV'  'AG003'
        'AG004_dH'  'AG004_dV'  'AG004'
        'AG005_dH'  'AG005_dV'  'AG005'
        'AG006_dH'  'AG006_dV'  'AG006'
        'AG007_dH'  'AG007_dV'  'AG007'
        'AG008_dH'  'AG008_dV'  'AG008'
        'AG009_dH'  'AG009_dV'  'AG009'
        'AG010_dH'  'AG010_dV'  'AG010'
        'AG011_dH'  'AG011_dV'  'AG011'
        'AG012_dH'  'AG012_dV'  'AG012'
        'AG013_dH'  'AG013_dV'  'AG013'
        'AG014_dH'  'AG014_dV'  'AG014'
        'AG015_dH'  'AG015_dV'  'AG015'
        'AG016_dH'  'AG016_dV'  'AG016'
        'AG017_dH'  'AG017_dV'  'AG017'
        'AG018_dH'  'AG018_dV'  'AG018'
        'AG019_dH'  'AG019_dV'  'AG019'
        'AG020_dH'  'AG020_dV'  'AG020'
        'AG021_dH'  'AG021_dV'  'AG021'
        'AG022_dH'  'AG022_dV'  'AG022'
        'AG023_dH'  'AG023_dV'  'AG023'
        'AG024_dH'  'AG024_dV'  'AG024'
        'AG025_dH'  'AG025_dV'  'AG025'
        'AG026_dH'  'AG026_dV'  'AG026'
        'AG027_dH'  'AG027_dV'  'AG027'
        'AG028_dH'  'AG028_dV'  'AG028'
        'AG029_dH'  'AG029_dV'  'AG029'
        'AG030_dH'  'AG030_dV'  'AG030'
        'AG031_dH'  'AG031_dV'  'AG031'
        'AG032_dH'  'AG032_dV'  'AG032'
        'AG065_dH'  'AG065_dV'  'AG065'
        'AG066_dH'  'AG066_dV'  'AG066'
        'AG067_dH'  'AG067_dV'  'AG067'
        'AG068_dH'  'AG068_dV'  'AG068'
        'AG069_dH'  'AG069_dV'  'AG069'
        'AG070_dH'  'AG070_dV'  'AG070'
        'AG071_dH'  'AG071_dV'  'AG071'
        'AG072_dH'  'AG072_dV'  'AG072'
        'AG073_dH'  'AG073_dV'  'AG073'
        'AG074_dH'  'AG074_dV'  'AG074'
        'AG075_dH'  'AG075_dV'  'AG075'
        'AG076_dH'  'AG076_dV'  'AG076'
        'AG077_dH'  'AG077_dV'  'AG077'
        'AG078_dH'  'AG078_dV'  'AG078'
        'AG079_dH'  'AG079_dV'  'AG079'
        'AG080_dH'  'AG080_dV'  'AG080'
        'AG081_dH'  'AG081_dV'  'AG081'
        'AG082_dH'  'AG082_dV'  'AG082'
        'AG083_dH'  'AG083_dV'  'AG083'
        'AG084_dH'  'AG084_dV'  'AG084'
        'AG085_dH'  'AG085_dV'  'AG085'
        'AG086_dH'  'AG086_dV'  'AG086'
        'AG087_dH'  'AG087_dV'  'AG087'
        'AG088_dH'  'AG088_dV'  'AG088'
        'AG089_dH'  'AG089_dV'  'AG089'
        'AG090_dH'  'AG090_dV'  'AG090'
        'AG091_dH'  'AG091_dV'  'AG091'
        'AG092_dH'  'AG092_dV'  'AG092'
        'AG093_dH'  'AG093_dV'  'AG093'
        'AG094_dH'  'AG094_dV'  'AG094'
        'AG095_dH'  'AG095_dV'  'AG095'
        'AG096_dH'  'AG096_dV'  'AG096'
        'AG122_dH'  'AG122_dV'  'AG122'
        'AG129_dH'  'AG129_dV'  'AG129'
        'AG130_dH'  'AG130_dV'  'AG130'
        'AG131_dH'  'AG131_dV'  'AG131'
        'AG132_dH'  'AG132_dV'  'AG132'
        'AG133_dH'  'AG133_dV'  'AG133'
        'AG134_dH'  'AG134_dV'  'AG134'
        'AG135_dH'  'AG135_dV'  'AG135'
        'AG136_dH'  'AG136_dV'  'AG136'
        'AG137_dH'  'AG137_dV'  'AG137'
        'AG138_dH'  'AG138_dV'  'AG138'
        'AG139_dH'  'AG139_dV'  'AG139'
        'AG140_dH'  'AG140_dV'  'AG140'
        'AG141_dH'  'AG141_dV'  'AG141'
        'AG142_dH'  'AG142_dV'  'AG142'
        'AG143_dH'  'AG143_dV'  'AG143'
        'AG144_dH'  'AG144_dV'  'AG144'
        'AG145_dH'  'AG145_dV'  'AG145'
        'AG146_dH'  'AG146_dV'  'AG146'
        'AG147_dH'  'AG147_dV'  'AG147'
        'AG148_dH'  'AG148_dV'  'AG148'
        'AG149_dH'  'AG149_dV'  'AG149'
        'AG150_dH'  'AG150_dV'  'AG150'
        'AG151_dH'  'AG151_dV'  'AG151'
        'AG152_dH'  'AG152_dV'  'AG152'
        'AG153_dH'  'AG153_dV'  'AG153'
        'AG154_dH'  'AG154_dV'  'AG154'
        'AG155_dH'  'AG155_dV'  'AG155'
        'AG156_dH'  'AG156_dV'  'AG156'
        'AG157_dH'  'AG157_dV'  'AG157'
        'AG158_dH'  'AG158_dV'  'AG158'
        'AG159_dH'  'AG159_dV'  'AG159'
        'AG160_dH'  'AG160_dV'  'AG160'
        'AG161_dH'  'AG161_dV'  'AG161'
        'AG162_dH'  'AG162_dV'  'AG162'
        'AG163_dH'  'AG163_dV'  'AG163'
        'AG164_dH'  'AG164_dV'  'AG164'
        'AG165_dH'  'AG165_dV'  'AG165'
        'AG193_dH'  'AG193_dV'  'AG193'
        'AG194_dH'  'AG194_dV'  'AG194'
        'AG195_dH'  'AG195_dV'  'AG195'
        'AG196_dH'  'AG196_dV'  'AG196'
        'AG197_dH'  'AG197_dV'  'AG197'
        'AG198_dH'  'AG198_dV'  'AG198'
        'AG199_dH'  'AG199_dV'  'AG199'
        'AG200_dH'  'AG200_dV'  'AG200'
        'AG201_dH'  'AG201_dV'  'AG201'
        'AG202_dH'  'AG202_dV'  'AG202'
        'AG203_dH'  'AG203_dV'  'AG203'
        'AG204_dH'  'AG204_dV'  'AG204'
        'AG205_dH'  'AG205_dV'  'AG205'
        'AG206_dH'  'AG206_dV'  'AG206'
        'AG207_dH'  'AG207_dV'  'AG207'
        'AG208_dH'  'AG208_dV'  'AG208'
        'AG209_dH'  'AG209_dV'  'AG209'
        'AG210_dH'  'AG210_dV'  'AG210'
        'AG211_dH'  'AG211_dV'  'AG211'
        'AG212_dH'  'AG212_dV'  'AG212'
        'AG213_dH'  'AG213_dV'  'AG213'
        'AG214_dH'  'AG214_dV'  'AG214'
        'AG215_dH'  'AG215_dV'  'AG215'
        'AG216_dH'  'AG216_dV'  'AG216'
        'AG217_dH'  'AG217_dV'  'AG217'
        'AG218_dH'  'AG218_dV'  'AG218'
        'AG219_dH'  'AG219_dV'  'AG219'
        'AG220_dH'  'AG220_dV'  'AG220'
        'AG221_dH'  'AG221_dV'  'AG221'
        'AG222_dH'  'AG222_dV'  'AG222'
        'AG223_dH'  'AG223_dV'  'AG223'
        'AG224_dH'  'AG224_dV'  'AG224'
        'AG225_dH'  'AG225_dV'  'AG225'
        'AG226_dH'  'AG226_dV'  'AG226'
        'AG227_dH'  'AG227_dV'  'AG227'
        'AG257_dH'  'AG257_dV'  'AG257'
        'AG258_dH'  'AG258_dV'  'AG258'
        'AG259_dH'  'AG259_dV'  'AG259'
        'AG260_dH'  'AG260_dV'  'AG260'
        'AG261_dH'  'AG261_dV'  'AG261'
        'AG262_dH'  'AG262_dV'  'AG262'
        'AG263_dH'  'AG263_dV'  'AG263'
        'AG264_dH'  'AG264_dV'  'AG264'
        'AG265_dH'  'AG265_dV'  'AG265'
        'AG266_dH'  'AG266_dV'  'AG266'
        'AG267_dH'  'AG267_dV'  'AG267'
        'AG268_dH'  'AG268_dV'  'AG268'
        'AG269_dH'  'AG269_dV'  'AG269'
        'AG270_dH'  'AG270_dV'  'AG270'
        'AG271_dH'  'AG271_dV'  'AG271'
        'AG272_dH'  'AG272_dV'  'AG272'
        'AG273_dH'  'AG273_dV'  'AG273'
        'AG274_dH'  'AG274_dV'  'AG274'
        'AG275_dH'  'AG275_dV'  'AG275'
        'AG276_dH'  'AG276_dV'  'AG276'
        'AG277_dH'  'AG277_dV'  'AG277'
        'AG278_dH'  'AG278_dV'  'AG278'
        'AG279_dH'  'AG279_dV'  'AG279'
        'AG280_dH'  'AG280_dV'  'AG280'
        'AG281_dH'  'AG281_dV'  'AG281'
        'AG282_dH'  'AG282_dV'  'AG282'
        'AG283_dH'  'AG283_dV'  'AG283'
        'AG284_dH'  'AG284_dV'  'AG284'
        'AG285_dH'  'AG285_dV'  'AG285'
        'AG286_dH'  'AG286_dV'  'AG286'
        'AG287_dH'  'AG287_dV'  'AG287'
        'AG288_dH'  'AG288_dV'  'AG288'
        'AG321_dH'  'AG321_dV'  'AG321'
        'AG322_dH'  'AG322_dV'  'AG322'
        'AG323_dH'  'AG323_dV'  'AG323'
        'AG324_dH'  'AG324_dV'  'AG324'
        'AG325_dH'  'AG325_dV'  'AG325'
        'AG326_dH'  'AG326_dV'  'AG326'
        'AG327_dH'  'AG327_dV'  'AG327'
        'AG328_dH'  'AG328_dV'  'AG328'
        'AG329_dH'  'AG329_dV'  'AG329'
        'AG330_dH'  'AG330_dV'  'AG330'
        'AG331_dH'  'AG331_dV'  'AG331'
        'AG332_dH'  'AG332_dV'  'AG332'
        'AG333_dH'  'AG333_dV'  'AG333'
        'AG334_dH'  'AG334_dV'  'AG334'
        'AG335_dH'  'AG335_dV'  'AG335'
        'AG336_dH'  'AG336_dV'  'AG336'
        'AG337_dH'  'AG337_dV'  'AG337'
        'AG338_dH'  'AG338_dV'  'AG338'
        'AG339_dH'  'AG339_dV'  'AG339'
        'AG340_dH'  'AG340_dV'  'AG340'
        'AG341_dH'  'AG341_dV'  'AG341'
        'AG342_dH'  'AG342_dV'  'AG342'
        'AG343_dH'  'AG343_dV'  'AG343'
        'AG344_dH'  'AG344_dV'  'AG344'
        'AG345_dH'  'AG345_dV'  'AG345'
        'AG346_dH'  'AG346_dV'  'AG346'
        'AG347_dH'  'AG347_dV'  'AG347'
        'AG348_dH'  'AG348_dV'  'AG348'
        'AG349_dH'  'AG349_dV'  'AG349'
        'AG350_dH'  'AG350_dV'  'AG350'
        'AG351_dH'  'AG351_dV'  'AG351'
        'AG352_dH'  'AG352_dV'  'AG352'
        'AG378_dH'  'AG378_dV'  'AG378'
        'AG385_dH'  'AG385_dV'  'AG385'
        'AG386_dH'  'AG386_dV'  'AG386'
        'AG387_dH'  'AG387_dV'  'AG387'
        'AG388_dH'  'AG388_dV'  'AG388'
        'AG389_dH'  'AG389_dV'  'AG389'
        'AG390_dH'  'AG390_dV'  'AG390'
        'AG391_dH'  'AG391_dV'  'AG391'
        'AG392_dH'  'AG392_dV'  'AG392'
        };
      yokogawa440_planar_combined = label(:,3);
      label = label(:,1:2);

    case {'eeg' 'electrode'}
      % there is no default set of electrode labels for all possible EEG systems
      % but nevertheless the requested input type should not result in an error
      label = {};

    otherwise
      error('the requested sensor type "%s" is not supported', type);

  end % switch

  % remember this set of labels to speed up subsequent function calls
  eval(sprintf('%s = label;', type));
  clear label

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch output
  case 'normal'
    % return labels as 2*Nx1 cell-array for planar systems or 3*Nx1 for neuromag306
    % return labels as   Nx1 cell-array for non-planar systems
    label = eval(type);

  case 'planarcombined'
    % return labels as Nx3 cell-array for the planar channels, 3rd column contains the combination
    planar    = eval(type);
    combined  = eval([type '_combined']);
    label     = [planar(:,1:2) combined]; % magnetometers are in the 3rd column for neuromag306

  otherwise
    error('unsupported output "%s"', output);

end

function vol = ea_ft_headmodel_simbio(geom, varargin)

% FT_HEADMODEL_SIMBIO creates a volume conduction model of the head
% using the finite element method (FEM) for EEG. This function takes
% as input a volumetric mesh (hexahedral or tetrahedral) and
% returns as output a volume conduction model which can be used to
% compute leadfields.
%
% This implements
%       ...
%
% Use as
%   vol = ft_headmodel_simbio(geom,'conductivity', conductivities, ...)
%
% The geom is given as a volumetric mesh, using ft_datatype_parcellation
%   geom.pos = vertex positions
%   geom.tet/geom.hex = list of volume elements
%   geom.tissue = tissue assignment for elements
%   geom.tissuelabel = labels correspondig to tissues
%
% Required input arguments should be specified in key-value pairs and have
% to include
%   conductivity   = vector containing tissue conductivities using ordered
%                    corresponding to geom.tissuelabel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To run this on Windows the following packages are necessary:
%
% Microsoft Visual C++ 2008 Redistributable
%
% Intel Visual Fortran Redistributables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% $Id: ft_headmodel_simbio.m 8445 2013-09-03 10:01:42Z johvor $


% get the optional arguments
conductivity    = ea_ft_getopt(varargin, 'conductivity');

% start with an empty volume conductor
geom = ea_ft_datatype_parcellation(geom);
vol = [];
if isfield(geom,'pos')
  vol.pos = geom.pos;
else
  error('Vertex field is required!')
end

if isfield(geom,'tet')
  vol.tet = geom.tet;
elseif isfield(geom,'hex')
  vol.hex = geom.hex;
else
  error('Connectivity information is required!')
end

if isfield(geom,'tissue')
  vol.tissue = geom.tissue;
else
  error('No element indices declared!')
end

if isempty(conductivity)
  error('No conductivity information!')
end

if length(conductivity) >= length(unique(vol.tissue))
  vol.cond = conductivity;
else
  error('Wrong conductivity information!')
end

if ~isfield(geom,'tissuelabel')
  numlabels = size(unique(geom.tissue),1);
  vol.tissuelabel = {};
  ulabel = unique(labels);
  for i = 1:numlabels
    vol.tissuelabel{i} = num2str(ulabel(i));
  end
else
  vol.tissuelabel = geom.tissuelabel;
end

vol.stiff = ea_sb_calc_stiff(vol);
vol.type = 'simbio';

function parcellation = ea_ft_datatype_parcellation(parcellation, varargin)

% FT_DATATYPE_PARCELLATION describes the FieldTrip MATLAB structure for parcellated
% cortex-based data and atlases. A parcellation can either be indexed or probabilistic
% (see below).
%
% A parcellation describes the tissue types for each of the surface elements.
% Parcellations are often, but not always labeled. A parcellatoin can be used to
% estimate the activity from MEG data in a known region of interest. A surface-based
% atlas is basically a very detailled parcellation with an anatomical label for each
% vertex.
%
% An example of a surface based Brodmann parcellation looks like this
%
%              pos: [8192x3]         positions of the vertices forming the cortical sheet
%              tri: [16382x3]        triangles of the cortical sheet
%         coordsys: 'ctf'            the (head) coordinate system in which the vertex positions are expressed
%             unit: 'mm'             the units in which the coordinate system is expressed
%         brodmann: [8192x1 uint8]   values from 1 to N, the value 0 means unknown
%    brodmannlabel: {Nx1 cell}
%
% An alternative representation of this parcellation is
%
%              pos: [8192x3]           positions of the vertices forming the cortical sheet
%              tri: [16382x3]          triangles of the cortical sheet
%         coordsys: 'ctf'              the (head) coordinate system in which the vertex positions are expressed
%             unit: 'mm'               the units in which the coordinate system is expressed
%  Brodmann_Area_1: [8192x1 logical]   binary map representing the voxels belonging to the specific area
%  Brodmann_Area_2: [8192x1 logical]   binary map representing the voxels belonging to the specific area
%  Brodmann_Area_3: [8192x1 logical]   binary map representing the voxels belonging to the specific area
%  ...
%
% The examples above demonstrate that a parcellation can be indexed, i.e. consisting of
% subsequent integer numbers (1, 2, ...) or probabilistic, consisting of real numbers
% ranging from 0 to 1 that represent probabilities between 0% and 100%. An extreme case
% is one where the probability is either 0 or 1, in which case the probability can be
% represented as a binary or logical array.
%
% The only difference to the source data structure is that the parcellation structure
% contains the additional fields xxx and xxxlabel. See FT_DATATYPE_SOURCE for further
% details.
%
% Required fields:
%   - pos
%
% Optional fields:
%   - tri, coordsys, unit
%
% Deprecated fields:
%   - none
%
% Obsoleted fields:
%   - none
%
% Revision history:
% (2012/latest) The initial version was defined in accordance with the representation of
% a voxel-based segmentation.
%
% See also FT_DATATYPE, FT_DATATYPE_SOURCE, FT_DATATYPE_SEGMENTATION

% Copyright (C) 2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_datatype_parcellation.m 10213 2015-02-11 19:38:33Z roboos $

% get the optional input arguments, which should be specified as key-value pairs
version           = ea_ft_getopt(varargin, 'version', 'latest');
parcellationstyle = ea_ft_getopt(varargin, 'parcellationstyle');  % can be indexed or probabilistic

if strcmp(version, 'latest')
  parcelversion = '2012';
  sourceversion = 'latest';
  clear version
else
  parcelversion = version;
  sourceversion = version;
  clear version
end

if isempty(parcellation)
  return;
end

switch parcelversion
  case '2012'

    if isfield(parcellation, 'pnt')
      parcellation.pos = parcellation.pnt;
      parcellation = rmfield(parcellation, 'pnt');
    end

    % convert the inside/outside fields, they should be logical rather than an index
    if isfield(parcellation, 'inside')
      parcellation = ea_fixinside(parcellation, 'logical');
    end

    dim = size(parcellation.pos,1);

    % make a list of fields that represent a parcellation
    fn = fieldnames(parcellation);
    fn = setdiff(fn, 'inside'); % exclude the inside field from any conversions
    sel = false(size(fn));
    for i=1:numel(fn)
      sel(i) = isnumeric(parcellation.(fn{i})) && numel(parcellation.(fn{i}))==dim;
    end
    % only consider numeric fields of the correct size
    fn = fn(sel);

    % determine whether the style of the input fields is probabilistic or indexed
    [indexed, probabilistic] = ea_determine_segmentationstyle(parcellation, fn, dim);

    % ignore the fields that do not contain a parcellation
    sel = indexed | probabilistic;
    fn            = fn(sel);
    indexed       = indexed(sel);
    probabilistic = probabilistic(sel);

    if ~any(probabilistic) && ~any(indexed)
      % rather than being described with a tissue label for each vertex
      % it can also be described with a tissue label for each surface or volme element
      for i = 1:length(fn)
        fname = fn{i};
        switch fname
          case 'tri'
            dim = size(parcellation.tri,1);
          case 'hex'
            dim = size(parcellation.hex,1);
          case 'tet'
            dim = size(parcellation.tet,1);
        end
      end
      [indexed, probabilistic] = ea_determine_segmentationstyle(parcellation, fn, dim);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ensure that the parcellation is internally consistent
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if any(probabilistic)
      parcellation = ea_fixsegmentation(parcellation, fn(probabilistic), 'probabilistic');
    end

    if any(indexed)
      parcellation = ea_fixsegmentation(parcellation, fn(indexed), 'indexed');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert the parcellation to the desired style
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmp(parcellationstyle, 'indexed') && any(probabilistic)
      parcellation  = convert_segmentationstyle(parcellation, fn(probabilistic), [dim 1], 'indexed');
    elseif strcmp(parcellationstyle, 'probabilistic') && any(indexed)
      parcellation  = convert_segmentationstyle(parcellation, fn(indexed), [dim 1], 'probabilistic');
    end % converting converting to desired style

  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for parcellation datatype', parcelversion);
end

% the parcellation is a speciat type of volume structure, so ensure that it also fulfills the requirements for that
parcellation = ea_ft_datatype_source(parcellation, 'version', sourceversion);

function source = ea_ft_datatype_source(source, varargin)

% FT_DATATYPE_SOURCE describes the FieldTrip MATLAB structure for data that is
% represented at the source level. This is typically obtained with a beamformer of
% minimum-norm source reconstruction using FT_SOURCEANALYSIS.
%
% An example of a source structure obtained after performing DICS (a frequency
% domain beamformer scanning method) is shown here
%
%           pos: [6732x3 double]       positions at which the source activity could have been estimated
%        inside: [6732x1 logical]      boolean vector that indicates at which positions the source activity was estimated
%           dim: [xdim ydim zdim]      if the positions can be described as a 3D regular grid, this contains the
%                                       dimensionality of the 3D volume
%     cumtapcnt: [120x1 double]        information about the number of tapers per original trial
%          time: 0.100                 the latency at which the activity is estimated (in seconds)
%          freq: 30                    the frequency at which the activity is estimated (in Hz)
%           pow: [6732x120 double]     the estimated power at each source position
%     powdimord: 'pos_rpt'             defines how the numeric data has to be interpreted,
%                                       in this case 6732 dipole positions x 120 repetitions (i.e. trials)
%           cfg: [1x1 struct]          the configuration used by the function that generated this data structure
%
% Required fields:
%   - pos
%
% Optional fields:
%   - time, freq, pow, coh, eta, mom, ori, cumtapcnt, dim, transform, inside, cfg, dimord, other fields with a dimord
%
% Deprecated fields:
%   - method, outside
%
% Obsoleted fields:
%   - xgrid, ygrid, zgrid, transform, latency, frequency
%
% Historical fields:
%   - avg, cfg, cumtapcnt, df, dim, freq, frequency, inside, method,
%   outside, pos, time, trial, vol, see bug2513
%
% Revision history:
%
% (2014) The subfields in the avg and trial fields are now present in the
% main structure, e.g. source.avg.pow is now source.pow. Furthermore, the
% inside is always represented as logical vector.
%
% (2011) The source representation should always be irregular, i.e. not
% a 3-D volume, contain a "pos" field and not contain a "transform".
%
% (2010) The source structure should contain a general "dimord" or specific
% dimords for each of the fields. The source reconstruction in the avg and
% trial substructures has been moved to the toplevel.
%
% (2007) The xgrid/ygrid/zgrid fields have been removed, because they are
% redundant.
%
% (2003) The initial version was defined
%
% See also FT_DATATYPE, FT_DATATYPE_COMP, FT_DATATYPE_DIP, FT_DATATYPE_FREQ,
% FT_DATATYPE_MVAR, FT_DATATYPE_RAW, FT_DATATYPE_SOURCE, FT_DATATYPE_SPIKE,
% FT_DATATYPE_TIMELOCK, FT_DATATYPE_VOLUME

% Copyright (C) 2013-2014, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_datatype_source.m 10265 2015-03-04 12:20:41Z jansch $

% FIXME: I am not sure whether the removal of the xgrid/ygrid/zgrid fields
% was really in 2007

% get the optional input arguments, which should be specified as key-value pairs
version = ea_ft_getopt(varargin, 'version', 'latest');

if strcmp(version, 'latest') || strcmp(version, 'upcoming')
  version = '2014';
end

if isempty(source)
  return;
end

% old data structures may use latency/frequency instead of time/freq. It is
% unclear when these were introduced and removed again, but they were never
% used by any fieldtrip function itself
if isfield(source, 'frequency'),
  source.freq = source.frequency;
  source      = rmfield(source, 'frequency');
end
if isfield(source, 'latency'),
  source.time = source.latency;
  source      = rmfield(source, 'latency');
end

switch version
  case '2014'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ensure that it has individual source positions
    source = ea_fixpos(source);

    % ensure that it is always logical
    source = ea_fixinside(source, 'logical');

    % remove obsolete fields
    if isfield(source, 'method')
      source = rmfield(source, 'method');
    end
    if isfield(source, 'transform')
      source = rmfield(source, 'transform');
    end
    if isfield(source, 'xgrid')
      source = rmfield(source, 'xgrid');
    end
    if isfield(source, 'ygrid')
      source = rmfield(source, 'ygrid');
    end
    if isfield(source, 'zgrid')
      source = rmfield(source, 'zgrid');
    end

    if isfield(source, 'avg') && isstruct(source.avg)
      % move the average fields to the main structure
      fn = fieldnames(source.avg);
      for i=1:length(fn)
        dat = source.avg.(fn{i});
        if isequal(size(dat), [1 size(source.pos,1)])
          source.(fn{i}) = dat';
        else
          source.(fn{i}) = dat;
        end
        clear dat
      end % j
      source = rmfield(source, 'avg');
    end

    if isfield(source, 'inside')
      % the inside is by definition logically indexed
      probe = find(source.inside, 1, 'first');
    else
      % just take the first source position
      probe = 1;
    end

    if isfield(source, 'trial') && isstruct(source.trial)
      npos = size(source.pos,1);

      % concatenate the fields for each trial and move them to the main structure
      fn = fieldnames(source.trial);

      for i=1:length(fn)
        % some fields are descriptive and hence identical over trials
        if strcmp(fn{i}, 'csdlabel')
          source.csdlabel = dat;
          continue
        end

        % start with the first trial
        dat    = source.trial(1).(fn{i});
        datsiz = ea_getdimsiz(source, fn{i});
        nrpt   = datsiz(1);
        datsiz = datsiz(2:end);


        if iscell(dat)
          datsiz(1) = nrpt; % swap the size of pos with the size of rpt
          val  = cell(npos,1);
          indx = find(source.inside);
          for k=1:length(indx)
            val{indx(k)}          = nan(datsiz);
            val{indx(k)}(1,:,:,:) = dat{indx(k)};
          end
          % concatenate all data as {pos}_rpt_etc
          for j=2:nrpt
            dat = source.trial(j).(fn{i});
            for k=1:length(indx)
              val{indx(k)}(j,:,:,:) = dat{indx(k)};
            end

          end % for all trials
          source.(fn{i}) = val;

        else
          % concatenate all data as pos_rpt_etc
          val = nan([datsiz(1) nrpt datsiz(2:end)]);
          val(:,1,:,:,:) = dat(:,:,:,:);
          for j=2:length(source.trial)
            dat = source.trial(j).(fn{i});
            val(:,j,:,:,:) = dat(:,:,:,:);
          end % for all trials
          source.(fn{i}) = val;

%         else
%           siz = size(dat);
%           if prod(siz)==npos
%             siz = [npos nrpt];
%           elseif siz(1)==npos
%             siz = [npos nrpt siz(2:end)];
%           end
%           val = nan(siz);
%           % concatenate all data as pos_rpt_etc
%           val(:,1,:,:,:) = dat(:);
%           for j=2:length(source.trial)
%             dat = source.trial(j).(fn{i});
%             val(:,j,:,:,:) = dat(:);
%           end % for all trials
%           source.(fn{i}) = val;

        end
      end % for each field

      source = rmfield(source, 'trial');

    end % if trial

    % ensure that it has a dimord (or multiple for the different fields)
    source = ea_fixdimord(source);


  case '2011'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ensure that it has individual source positions
    source = ea_fixpos(source);

    % remove obsolete fields
    if isfield(source, 'xgrid')
      source = rmfield(source, 'xgrid');
    end
    if isfield(source, 'ygrid')
      source = rmfield(source, 'ygrid');
    end
    if isfield(source, 'zgrid')
      source = rmfield(source, 'zgrid');
    end
    if isfield(source, 'transform')
      source = rmfield(source, 'transform');
    end

    % ensure that it has a dimord (or multiple for the different fields)
    source = ea_fixdimord(source);

  case '2010'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ensure that it has individual source positions
    source = ea_fixpos(source);

    % remove obsolete fields
    if isfield(source, 'xgrid')
      source = rmfield(source, 'xgrid');
    end
    if isfield(source, 'ygrid')
      source = rmfield(source, 'ygrid');
    end
    if isfield(source, 'zgrid')
      source = rmfield(source, 'zgrid');
    end

    % ensure that it has a dimord (or multiple for the different fields)
    source = ea_fixdimord(source);

  case '2007'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ensure that it has individual source positions
    source = ea_fixpos(source);

    % remove obsolete fields
    if isfield(source, 'dimord')
      source = rmfield(source, 'dimord');
    end
    if isfield(source, 'xgrid')
      source = rmfield(source, 'xgrid');
    end
    if isfield(source, 'ygrid')
      source = rmfield(source, 'ygrid');
    end
    if isfield(source, 'zgrid')
      source = rmfield(source, 'zgrid');
    end

  case '2003'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(source, 'dimord')
      source = rmfield(source, 'dimord');
    end

    if ~isfield(source, 'xgrid') || ~isfield(source, 'ygrid') || ~isfield(source, 'zgrid')
      if isfield(source, 'dim')
        minx = min(source.pos(:,1));
        maxx = max(source.pos(:,1));
        miny = min(source.pos(:,2));
        maxy = max(source.pos(:,2));
        minz = min(source.pos(:,3));
        maxz = max(source.pos(:,3));
        source.xgrid = linspace(minx, maxx, source.dim(1));
        source.ygrid = linspace(miny, maxy, source.dim(2));
        source.zgrid = linspace(minz, maxz, source.dim(3));
      end
    end

  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for source datatype', version);
end

function [source] = ea_fixinside(source, opt)

% FIXINSIDE ensures that the region of interest (which is indicated by the
% field "inside") is consistently defined for source structures and volume
% structures. Furthermore, it solves backward compatibility problems.
%
% Use as
%   [source] = fixinside(source, 'logical');
% or
%   [source] = fixinside(source, 'index');

% Copyright (C) 2006, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: fixinside.m 9664 2014-06-22 07:06:29Z roboos $


if nargin<2
  opt = 'logical';
end

if ~isfield(source, 'inside')
  if isfield(source, 'pos')
    % assume that all positions are inside the region of interest
    source.inside  = [1:size(source.pos,1)]';
    source.outside = [];
  elseif isfield(source, 'dim')
    source.inside  = [1:prod(source.dim)]';
    source.outside = [];
  end
end

if ~isfield(source, 'inside')
  % nothing to do
  return;
end

% determine the format
if isa(source.inside, 'logical')
  logicalfmt = 1;
elseif all(source.inside(:)==0 | source.inside(:)==1)
  source.inside = logical(source.inside);
  logicalfmt = 1;
else
  logicalfmt = 0;
end

if ~logicalfmt && strcmp(opt, 'logical')
  % convert to a logical array
  if ~isfield(source, 'outside')
    source.outside = [];
  end
  inside(source.inside)  = (1==1);  % true
  inside(source.outside) = (1==0);  % false
  source.inside = inside(:);
  if isfield(source, 'outside')
    source = rmfield(source, 'outside');
  end
elseif logicalfmt && strcmp(opt, 'index')
  % convert to a vectors with indices
  tmp = source.inside;
  source.inside  = find( tmp(:));
  source.outside = find(~tmp(:));
else
  % nothing to do
end

function dimsiz = ea_getdimsiz(data, field)

% GETDIMSIZ
%
% Use as
%   dimsiz = getdimsiz(data, field)
%
% See also GETDIMORD

if ~isfield(data, field) && isfield(data, 'avg') && isfield(data.avg, field)
  field = ['avg.' field];
elseif ~isfield(data, field) && isfield(data, 'trial') && isfield(data.trial, field)
  field = ['trial.' field];
elseif ~isfield(data, field)
  error('field "%s" not present in data', field);
end

if strncmp(field, 'avg.', 4)
  prefix = [];
  field = field(5:end); % strip the avg
  data.(field) = data.avg.(field); % move the avg into the main structure
  data = rmfield(data, 'avg');
elseif strncmp(field, 'trial.', 6)
  prefix = numel(data.trial);
  field = field(7:end); % strip the trial
  data.(field) = data.trial(1).(field); % move the first trial into the main structure
  data = rmfield(data, 'trial');
else
  prefix = [];
end

dimsiz = ea_cellmatsize(data.(field));

% add nrpt in case of source.trial
dimsiz = [prefix dimsiz];

function siz = ea_cellmatsize(x)
if iscell(x)
  cellsize = numel(x);          % the number of elements in the cell-array
  [dum, indx] = max(cellfun(@numel, x));
  matsize = size(x{indx});      % the size of the content of the cell-array
  siz = [cellsize matsize];     % concatenate the two
else
  siz = size(x);
end

function [data] = ea_fixdimord(data)

% FIXDIMORD ensures consistency between the dimord string and the axes
% that describe the data dimensions. The main purpose of this function
% is to ensure backward compatibility of all functions with data that has
% been processed by older FieldTrip versions
%
% Use as
%   [data] = fixdimord(data)
% This will modify the data.dimord field to ensure consistency.
% The name of the axis is the same as the name of the dimord, i.e. if
% dimord='freq_time', then data.freq and data.time should be present.
%
% The default dimensions in the data are described by
%  'time'
%  'freq'
%  'chan'
%  'chancmb'
%  'refchan'
%  'subj'
%  'rpt'
%  'rpttap'
%  'pos'
%  'ori'
%  'rgb'
%  'comp'
%  'voxel'

% Copyright (C) 2009-2014, Robert Oostenveld, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: fixdimord.m 9972 2014-11-19 08:09:34Z roboos $

% if nargin<2, keepsourcedimord = 0; end
%
% if any(ft_datatype(data, {'source', 'volume'})) && isfield(data, 'dimord') && ~keepsourcedimord
%   % the old source data representation does not have a dimord, whereas the new source data representation does have a dimord
%   warning(sprintf('removing dimord "%s" from source representation data', data.dimord));
%   data = rmfield(data, 'dimord');
%   return
% else
%   % it is ok
%   return
% end

if ~isfield(data, 'dimord')
  if ft_datatype(data, 'raw')
    % it is raw data, which does not have a dimord -> this is ok
    return
  elseif ft_datatype(data, 'comp')
    % it is component data, which resembles raw data -> this is ok
    return
  elseif ft_datatype(data, 'volume')
    % it is volume data, which does not have a dimord -> this is ok
    return
  else
    fn = fieldnames(data);
    sel = true(size(fn));
    for i=1:length(fn)
      sel(i) = contains(fn{i}, 'dimord');
    end
    df = fn(sel);

    if isempty(df)
      if ft_datatype(data, 'source') || ft_datatype(data, 'parcellation')
        % it is old-style source data -> this is ok
        % ft_checkdata will convert it to new-style
        return
      else
        error('the data does not contain a dimord, but it also does not resemble raw or component data');
      end
    end

    % use this function recursively on the XXXdimord fields
    for i=1:length(df)
      data.dimord = data.(df{i});
      data = fixdimord(data);
      data.(df{i}) = data.dimord;
      data = rmfield(data, 'dimord');
    end
    % after the recursive call it should be ok
    return
  end
end

if strcmp(data.dimord, 'voxel')
  % this means that it is position
  data.dimord = 'pos';
end

dimtok = tokenize(data.dimord, '_');
if strncmp('{pos_pos}', data.dimord, 9)
  % keep these together for bivariate source structures
  dimtok = {'{pos_pos}', dimtok{3:end}};
end

for i=1:length(dimtok)
  switch dimtok{i}
    case {'tim' 'time' 'toi' 'latency'}
      dimtok{i} = 'time';

    case {'frq' 'freq' 'foi' 'frequency'}
      dimtok{i} = 'freq';

    case {'sgn' 'label' 'chan'}
      dimtok{i} = 'chan';

    case {'rpt' 'trial'}
      dimtok{i} = 'rpt';

    case {'subj' 'subject'}
      dimtok{i} = 'subj';

    case {'comp'}
      % don't change, it is ok

    case {'sgncmb' 'labelcmb' 'chancmb'}
      dimtok{i} = 'chancmb';

    case {'rpttap'}
      % this is a 2-D field, coding trials and tapers along the same dimension
      % don't change, it is ok

    case {'refchan'}
      % don't change, it is ok

    case {'ori'}
      % don't change, it is ok

    case {'rgb'}
      % don't change, it is ok

    case {'voxel' 'vox' 'repl' 'wcond'}
      % these are used in some fieldtrip functions, but are not considered standard
      warning_once('unexpected dimord "%s"', data.dimord);

    case {'pos'}
      % this is for source data on a 3-d grid, a cortical sheet, or unstructured positions

    case {'{pos}' '{pos}_rpt' '{pos}_rpttap'}
      % this is for source data on a 3-d grid, a cortical sheet, or unstructured positions
      % the data itself is represented in a cell-array, e.g. source.mom or source.leadfield

    case {'{pos_pos}'}
      % this is for bivariate source data on a 3-d grid, a cortical sheet, or unstructured positions

    otherwise
      error(sprintf('unexpected dimord "%s"', data.dimord));

  end % switch dimtok
end % for length dimtok

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(data, 'tim'),         data.time      = data.tim         ; data = rmfield(data, 'tim')        ; end
if isfield(data, 'toi'),         data.time      = data.toi         ; data = rmfield(data, 'toi')        ; end
if isfield(data, 'latency'),     data.time      = data.latency     ; data = rmfield(data, 'latency')    ; end
if isfield(data, 'frq'),         data.freq      = data.frq         ; data = rmfield(data, 'frq')        ; end
if isfield(data, 'foi'),         data.freq      = data.foi         ; data = rmfield(data, 'foi')        ; end
if isfield(data, 'frequency'),   data.freq      = data.frequency   ; data = rmfield(data, 'frequency')  ; end
if isfield(data, 'sgn'),         data.label     = data.sgn         ; data = rmfield(data, 'sgn')        ; end
if isfield(data, 'chan'),        data.label     = data.chan        ; data = rmfield(data, 'chan')       ; end
% if isfield(data, 'trial'),         data.rpt     = data.trial         ; data = rmfield(data, 'trial')        ; end  % DO NOT CONVERT -> this is an exception
if isfield(data, 'subject'),     data.subj      = data.subject     ; data = rmfield(data, 'subject')    ; end
if isfield(data, 'sgncmb'),      data.labelcmb  = data.sgncmb      ; data = rmfield(data, 'sgncmb')     ; end
if isfield(data, 'chancmb'),     data.labelcmb  = data.chancmb     ; data = rmfield(data, 'chancmb')    ; end

% ensure that it is a column
if isfield(data, 'label')
  data.label = data.label(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if isfield(data, 'trial')
%   mat = data.trial;
% elseif isfield(data, 'individual')
%   mat = data.individual;
% elseif isfield(data, 'avg')
%   mat = data.avg;
% elseif isfield(data, 'crsspctrm')
%   mat = data.crsspctrm;
% elseif isfield(data, 'powspctrm')
%   mat = data.powspctrm;
% elseif isfield(data, 'fourierspctrm')
%   mat = data.fourierspctrm;
% end
%
% add the descriptive axis for each dimension
% for i=1:length(dimtok)
%   if isfield(data, dimtok{i})
%     % the dimension is already described with its own axis
%     % data = setfield(data, dimtok{i}, getfield(data, dimtok{i}));
%   else
%     % add an axis to the output data
%     data = setfield(data, dimtok{i}, 1:size(mat,i));
%   end
% end

% undo the tokenization
data.dimord = dimtok{1};
for i=2:length(dimtok)
  data.dimord = [data.dimord '_' dimtok{i}];
end

function source = ea_fixpos(source)
if ~isfield(source, 'pos')
  if isfield(source, 'xgrid') && isfield(source, 'ygrid') && isfield(source, 'zgrid')
    source.pos = ea_grid2pos(source.xgrid, source.ygrid, source.zgrid);
  elseif isfield(source, 'dim') && isfield(source, 'transform')
    source.pos = ea_dim2pos(source.dim, source.transform);
  else
    error('cannot reconstruct individual source positions');
  end
end

function pos = ea_grid2pos(xgrid, ygrid, zgrid)
[X, Y, Z] = ndgrid(xgrid, ygrid, zgrid);
pos = [X(:) Y(:) Z(:)];

function pos = ea_dim2pos(dim, transform)
[X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
pos = [X(:) Y(:) Z(:)];
pos = ea_ft_warp_apply(transform, pos, 'homogenous');

function [indexed, probabilistic] = ea_determine_segmentationstyle(segmentation, fn, dim)

% DETERMINE_SEGMENTATIONSTYLE is a helper function that determines the type of segmentation
% contained in each of the fields. It is used by ft_datatype_segmentation and
% ft_datatype_parcellation.
%
% See also FIXSEGMENTATION, CONVERT_SEGMENTATIONSTYLE

indexed       = false(size(fn));
probabilistic = false(size(fn));

% determine for each of the fields whether it is probabilistic, indexed or something else
for i=1:numel(fn)
  if numel(segmentation.(fn{i}))~=prod(dim)
    % this does not look like a segmentation
    continue
  elseif strcmp(fn{i}, 'anatomy')
    % this should not be interpreted as segmentation, also not when it is a uint8 or uint16 representation
    continue
  else
    if isfield(segmentation, [fn{i} 'label'])
      % the xxxlabel field exists, which only makes sense for an indexed representation
      probabilistic(i) = false;
      indexed(i)       = true;
    else
      % this looks like a segmentation
      tmp = segmentation.(fn{i});
      tmp = tmp(:);       % convert to vector
      sel = isnan(tmp);   % find NaN values
      if any(sel)
        % note that the the finding and removing of NaNs have been separated to speed up the code
        tmp = tmp(~sel);  % remove NaN values
      end
      clear sel
      probabilistic(i) =  islogical(tmp) || all(tmp>=-0.001 & tmp<=1.001); % allow some roundoff error
      indexed(i)       = ~islogical(tmp) && all(abs(tmp - round(tmp))<1000*eps);

      if probabilistic(i) && indexed(i)
        % the xxxlabel does not exist, so treat it as a probabilistic representation
        probabilistic(i) = true;
        indexed(i)       = false;
      end
    end
  end
end % for each of the fields

function rows = ea_sb_sparse_to_mat(diinsy)

% SB_SPARSE_TO_MAT
%
% $Id: sb_sparse_to_mat.m 8776 2013-11-14 09:04:48Z roboos $

rows = zeros(max(diinsy),1);
rows(diinsy) = 1;
rows = [1;rows];
rows = cumsum(rows);
rows = rows(1:end-1);

function [stiff, diinsy, cols, sysmat] = ea_sb_calc_stiff(vol)

% SB_CALC_STIFF
%
% $Id: sb_calc_stiff.m 8776 2013-11-14 09:04:48Z roboos $

if(~(size(vol.pos,2)==3))
    if(size(vol.pos,1)==3)
        node = vol.pos';
        warning('Dimension of vol.pos should be #nodes x 3!')
    else
        error('vol.pos has wrong dimension!')
    end
else
    node = vol.pos;
end
npnt = size(node,1);
npnt = int32(npnt);

if isfield(vol,'tet')
    if size(vol.tet,1) == 4
        mele = size(vol.tet,1);
        elem = vol.tet;
    elseif size(vol.tet,2) == 4
        mele = size(vol.tet,2);
        elem = vol.tet';
    else
        error('vol.tet has wrong dimensions!')
    end
    elem = [elem; zeros(4,size(elem,2))];
elseif isfield(vol,'hex')
    if size(vol.hex,1) == 8
        mele = size(vol.hex,1);
        elem = vol.hex;
    elseif size(vol.hex,2) == 8
        mele = size(vol.hex,2);
        elem = vol.hex';
    else
        error('vol.hex has wrong dimensions!')
    end
else
    error('Could not find connectivity information!')
end

if min(min(elem(1:mele,:))) == 0
    elem = elem + 1;
    warning('Numbering of nodes in vol.tet/vol.hex must start at 1 (Fortran numbering)!')
elseif min(min(elem(1:mele,:))) < 0
    error('No negative indices for conectivity information allowed!')
end

if isfield(vol,'cond') && isfield(vol,'tissue') && isfield(vol,'tissuelabel')
    if length(vol.tissuelabel) == length(vol.cond)
         if length(vol.tissue) == size(elem,2)
             cond = zeros(size(elem,2),1);
             numlabels = length(vol.tissuelabel);
             for i=1:numlabels
                 cond(vol.tissue == i) = vol.cond(i);
             end
        else
            error('Dimensions of vol.tet or vol.hex and vol.tissue do not fit!');
        end
    else
        error('Dimensions of vol.cond and entries of vol.tissuelabel do not fit!');
    end
end

mele = int32(mele);
elem = int32(elem);

% check whether the nodes have right orientation

if isfield(vol,'tet')
    if ~ea_sb_test_ori(node,elem(1:4,:)')
        error('Elements have wrong orientation, consider exchanging node 3 and 4');
        return;
    end
elseif isfield(vol,'hex')
    if ~ea_sb_test_ori(node,elem')
        error('Elements have wrong orientation or are degenerated');
        return
    end
end

try
    [diinsy,cols,sysmat] = ea_calc_stiff_matrix_val_wrapper(node,elem,cond,mele);
    ea_delete([pwd, filesep, 'fort.6']);
catch err
    if ispc && strcmp(err.identifier,'MATLAB:invalidMEXFile')
        error('Error executing mex-file. Microsoft Visual C++ 2008 Redistributables and Intel Visual Fortran Redistributables are required.')
    else
        rethrow(err)
    end
end
npnt = double(npnt);
diinsy = double(diinsy);
cols = double(cols);
rows = ea_sb_sparse_to_mat(diinsy);
stiff = sparse(rows,cols,sysmat,npnt,npnt,length(sysmat));

function err = ea_sb_test_ori(pnt,elem)
err = 1;
if(size(elem,2) == 4)
    det = sum(cross(pnt(elem(:,2),:)-pnt(elem(:,1),:),pnt(elem(:,4),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,3),:)-pnt(elem(:,1),:)),2);
    if length(find(det <= 0)) > 0
        err = 0;
    end
elseif(size(elem,2) == 8)
    det1 = sum(cross(pnt(elem(:,6),:)-pnt(elem(:,1),:),pnt(elem(:,8),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,5),:)-pnt(elem(:,1),:)),2);
    det2 = sum(cross(pnt(elem(:,3),:)-pnt(elem(:,1),:),pnt(elem(:,6),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,2),:)-pnt(elem(:,1),:)),2);
    det3 = sum(cross(pnt(elem(:,8),:)-pnt(elem(:,1),:),pnt(elem(:,3),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,4),:)-pnt(elem(:,1),:)),2);
    det4 = sum(cross(pnt(elem(:,8),:)-pnt(elem(:,3),:),pnt(elem(:,6),:)-pnt(elem(:,3),:),2).*(pnt(elem(:,7),:)-pnt(elem(:,3),:)),2);
    if (length(find(det1 <= 0)) > 0 || length(find(det2 <= 0)) > 0 || length(find(det3 <= 0)) > 0 || length(find(det4 <= 0)) > 0)
        err = 0;
    end
else
    error('Invalid number of nodes per element!');
end

