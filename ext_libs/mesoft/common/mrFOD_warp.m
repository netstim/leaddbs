function mr = mr_warp(idname,dname,mrname,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  mr = mr_FODwarp(idname,dname,ftrname)
%
% computes a nonrigid deformation of a FOD (mrstruct) according to the 
% deformation fields computed by SPM newsegment.
% The mrstruct is 4D [w h d  numorientations]
%
%
% PARAMETERS
%  idname - filename of inverse deformation (iy*.mat)
%  dname - filename of deformation (y*.mat)
%  mrname - filename of mrstruct or mrstruct itself
%  varargin - 'modulate', 0/1/2  do modulation by  1 or |Jn|^4/|J|^2 or |Jn|^3/|J| 
%           - 'outdirs', directions in output space ([3 n] array)
%           - 'method', interpolation method like in spm_sample_vol
%           - 'noutdir', number of output directions if outdirs is not defined
% OUTPUT
%  mr - deformed FOD as mrstruct
%
% Marco Reisert 2012, UKLFR
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


modulate = 1;
outdirs = [];
method = 1;
noutdir = 64;


if isstr(mrname),
    display('loading mrstruct');
    mr = mrstruct_read(mrname);
else,
    mr = mrname;
end;
ndir = size(mr.user.bDir,2);

try 
    sinterp = load(sprintf('sinterp%istruct.mat',ndir));
catch
    warning('sinterpstruct not existing, creating one');
    sinterp = sphereInterpolLUT(mr.user.bDir');
end;


for k = 1:2:length(varargin),
    if not(exist(varargin{k})),
        display(['invalid parameter: ' varargin{k}]);
    else,
        eval(sprintf('%s=varargin{k+1};',varargin{k}));
    end;
end;

if isempty(outdirs),
    sinterp_out = load(sprintf('sinterp%istruct.mat',noutdir));
    outdirs = sinterp_out.bDir;
end;

mr = mr_warp_int(dname,idname,mr,outdirs,sinterp,modulate,method);
mr.dataAy = mr.dataAy(:,:,:,1:noutdir);
mr.user.bDir = mr.user.bDir(:,1:noutdir);


return;


function mr_res= mr_warp_int(dname,idname,mr,outdirs,sinterp,modulate,method)

display('loading deformation fields');
invdef = load_def(idname);
sdef = load_def(dname);

bDir = outdirs;
sym = true;




display('warping dirs');
invjac = single(inv_jacobian(invdef,sdef,mr,method));
sz = [size(mr.dataAy,1) size(mr.dataAy,2) size(mr.dataAy,3)];
tdirfield = single(zeros([3,size(bDir,2),sz]));
for k = 1:3,
    for j = 1:3,        
        tdirfield(k,:,:,:,:)  =  squeeze(tdirfield(k,:,:,:,:)) + repmat( reshape(invjac(:,:,:,k,j),[1 sz]),[size(bDir,2) 1]  ).*squeeze(repmat(bDir(j,:) , [1 1 sz]));
    end;
end;   
dat = single(permute(mr.dataAy,[4 1 2 3]));  

if size(dat,1) == sinterp.numpoints,
    xxx = evalSinterp(tdirfield(:,:,:),dat(:,:),sinterp);
    sym = false;
elseif size(dat,1)*2 == sinterp.numpoints,
    xxx = evalSinterp(tdirfield(:,:,:),[dat(:,:) ; dat(:,:)],sinterp);
else
    display('something wrong with directions');
    return
end;

mr_res = mr; 
mr_res.user.bDir = bDir; 
mr_res.user.sym = sym;
mr_res.dataAy = permute(reshape(xxx,[size(bDir,2) sz]),[2 3 4 1]);

if modulate == 1, % surface integral preserving
    
    display('modulate with surface integral preserving factor');
    iJn4 = squeeze(single(sum(tdirfield.^2)).^2);
    detiJ2 = (invjac(:,:,:,1,1).*invjac(:,:,:,2,2).*invjac(:,:,:,3,3) + invjac(:,:,:,2,1).*invjac(:,:,:,3,2).*invjac(:,:,:,1,3) + invjac(:,:,:,3,1).*invjac(:,:,:,2,2).*invjac(:,:,:,2,3) + ...
          - invjac(:,:,:,1,3).*invjac(:,:,:,2,2).*invjac(:,:,:,3,1) - invjac(:,:,:,2,3).*invjac(:,:,:,3,2).*invjac(:,:,:,1,1) - invjac(:,:,:,3,3).*invjac(:,:,:,1,2).*invjac(:,:,:,2,1)).^2;

    for k = 1:size(mr_res.dataAy,4),
        mr_res.dataAy(:,:,:,k) = squeeze(mr_res.dataAy(:,:,:,k)) .* (squeeze(iJn4(k,:,:,:))./(eps+detiJ2));
    end;
 
elseif modulate == 2, % fiber density preserving
    
    display('modulate with fiber density preserving factor');
    iJn3 = abs(squeeze(single(sum(tdirfield.^2)))).^(3/2);
    detiJ = abs((invjac(:,:,:,1,1).*invjac(:,:,:,2,2).*invjac(:,:,:,3,3) + invjac(:,:,:,2,1).*invjac(:,:,:,3,2).*invjac(:,:,:,1,3) + invjac(:,:,:,3,1).*invjac(:,:,:,2,2).*invjac(:,:,:,2,3) + ...
          - invjac(:,:,:,1,3).*invjac(:,:,:,2,2).*invjac(:,:,:,3,1) - invjac(:,:,:,2,3).*invjac(:,:,:,3,2).*invjac(:,:,:,1,1) - invjac(:,:,:,3,3).*invjac(:,:,:,1,2).*invjac(:,:,:,2,1)));

    for k = 1:size(mr_res.dataAy,4),
        mr_res.dataAy(:,:,:,k) = squeeze(mr_res.dataAy(:,:,:,k)) .* (squeeze(iJn3(k,:,:,:))./(eps+detiJ));
    end;
end;

display('warping field');
mr_res.dataAy(isnan(mr_res.dataAy(:))) = 0;
mr_res.dataAy = single(forward(invdef,sdef,mr_res,method));


q = diag([-1 -1 1 1]);  % matrix for switching world coordinate conventions of nifti and mrstruct
nx = size(sdef.img,2);
p = [0 -1 0 nx-1; 1 0 0 0; 0 0 1 0; 0 0 0 1];   % matrix for switching voxel coordinate conventions of nifti and mrstruct
mr_res.edges = q*sdef.mat; %*inv(p);


mr_res.user.mFOD = squeeze(mean(mr_res.dataAy,4));
mr_res.user.mask = mr_res.user.mFOD > 0;

return



% forward warp by y-field of 4D-mrstruct

function warped = forward(invdef,sdef,mr,method)



q = diag([-1 -1 1 1]);  % matrix for switching world coordinate conventions of nifti and mrstruct

col = @(x) x(:);

def(1,:) = col(sdef.img(:,:,:,1));
def(2,:) = col(sdef.img(:,:,:,2));
def(3,:) = col(sdef.img(:,:,:,3));
T = inv(mr.edges)*q;
Tdef = double(T(1:3,1:3)*def + repmat(T(1:3,4),[1 size(def,2)]));
sz = size(sdef.img(:,:,:,1));

warped = zeros([size(sdef.img(:,:,:,1)) size(mr.dataAy,4)]);

k = 1;
warped(:,:,:,k) = reshape(spm_sample_vol(mr.dataAy(:,:,:,k),Tdef(1,:)+1,Tdef(2,:)+1,Tdef(3,:)+1,[method nan]),sz);
valid = (imfilter(double(not(isnan(warped(:,:,:,k)))),ones(5,5,5))==125);
warped(:,:,:,k) = warped(:,:,:,k).*valid;


for k = 2:size(mr.dataAy,4),
    warped(:,:,:,k) = reshape(spm_sample_vol(mr.dataAy(:,:,:,k),Tdef(1,:)+1,Tdef(2,:)+1,Tdef(3,:)+1,[method nan]),sz).*valid;
    fprintf('.');
end;    

fprintf('\n');
    








%%%%%%%%%% reslices iy-deformationfield into mrstruct-voxelspace

function Tidef3 = reslice_invdef(invdef,sdef,mr,method)
 
q = diag([-1 -1 1 1]);  % matrix for switching world coordinate conventions of nifti and mrstruct

A = inv(invdef.mat)*q*mr.edges;
[X Y Z] = ndgrid(0:size(mr.dataAy,1)-1,0:size(mr.dataAy,2)-1,0:size(mr.dataAy,3)-1);
C = [X(:)' ; Y(:)' ; Z(:)'];
AC = A(1:3,1:3)*C + repmat(A(1:3,4),[1 size(C,2)]);
idef(1,:) = spm_sample_vol(invdef.img(:,:,:,1),AC(1,:)+1,AC(2,:)+1,AC(3,:)+1,[method nan]);
idef(2,:) = spm_sample_vol(invdef.img(:,:,:,2),AC(1,:)+1,AC(2,:)+1,AC(3,:)+1,[method nan]);
idef(3,:) = spm_sample_vol(invdef.img(:,:,:,3),AC(1,:)+1,AC(2,:)+1,AC(3,:)+1,[method nan]);

W = inv(sdef.mat);
Tidef = double(W(1:3,1:3)*idef + repmat(W(1:3,4),[1 size(idef,2)]));

sz = size(mr.dataAy);
Tidef3(:,:,:,1) = reshape(Tidef(1,:),sz(1:3));
Tidef3(:,:,:,2) = reshape(Tidef(2,:),sz(1:3));
Tidef3(:,:,:,3) = reshape(Tidef(3,:),sz(1:3));


% backward warp by iy

function iwarped = backward(invdef,sdef,mr,img,method)

Tidef = reslice_invdef(invdef,sdef,mr);
iwarped = reshape(spm_sample_vol(img,Tidef(:,:,:,1)+1,Tidef(:,:,:,2)+1,Tidef(:,:,:,3)+1,[method nan]),size(mr.dataAy));


% compute jacobian by finite differences

function jac = inv_jacobian(invdef,sdef,mr,method)
Tidef = reslice_invdef(invdef,sdef,mr,method);

for k = 1:3,
    jac(:,:,:,1,k) =  (circshift(Tidef(:,:,:,k),[-1 0 0]) - circshift(Tidef(:,:,:,k),[1 0 0]))/2;
    jac(:,:,:,2,k) =  (circshift(Tidef(:,:,:,k),[0 -1 0]) - circshift(Tidef(:,:,:,k),[0 1 0]))/2;
    jac(:,:,:,3,k) =  (circshift(Tidef(:,:,:,k),[0 0 -1]) - circshift(Tidef(:,:,:,k),[0 0 1]))/2;
end;




% Load a deformation field saved as an image

function iy = load_def(name)

P      = [repmat(name,3,1), [',1,1';',1,2';',1,3']];
V      = spm_vol(P);
for k = 1:3,
    Def(:,:,:,k) = spm_load_float(V(k));
end;
mat    = V(1).mat;

iy.img = Def;

corrMy= diag(ones(1, 4)); corrMy(1:3, 4)= 1;  % switch to zero-based voxels
iy.mat = mat * corrMy;




return

