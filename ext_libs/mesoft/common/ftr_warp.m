function ftr = ftr_warp(idname,dname,ftrname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ftr = ftr_warp(idname,dname,ftrname)
%
% computes a nonrigid deformation of a _FTR.mat according to the 
% deformation fields computed by SPM newsegment
%
% PARAMETERS
%  idname - filename of inverse deformation (iy*.mat)
%  dname - filename of deformation (y*.mat)
%  ftrname - filename of FTRstruct
% OUTPUT
%  ftr - deformed FTRstruct
%
% Marco Reisert 2012, UKLFR
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



ftr = ftr_warp_int(dname,idname,ftrname);



function ftr = ftr_warp_int(dname,idname,ftrname)

display('loading deformation fields');
idef = load_def(idname);
sdef = load_def(dname);

display('loading fibers');
ftr = ftrstruct_read(ftrname);

display('deforming fibers');
q = diag([-1 -1 1 1]);  % matrix for switching world coordinate conventions of nifti and mrstruct
fibs = ftr.curveSegCell;
t = inv(idef.mat)*q*ftr.hMatrix;     % from mrs-voxelcoordinates of FTR to nifti-voxelcoordinates of deformationfield
tfibs = cellfun(@(x) (x-1)*t(1:3,1:3)' + repmat(t(1:3,4)',[size(x,1) 1]) ,fibs,'uniformoutput',false);

% applying deformation
fibo = double(cat(1,tfibs{:})); % to do it all at once with spm_sample_vol
nfibx = spm_sample_vol(idef.img(:,:,:,1),fibo(:,1),fibo(:,2),fibo(:,3),[1 nan]);
nfiby = spm_sample_vol(idef.img(:,:,:,2),fibo(:,1),fibo(:,2),fibo(:,3),[1 nan]);
nfibz = spm_sample_vol(idef.img(:,:,:,3),fibo(:,1),fibo(:,2),fibo(:,3),[1 nan]);
nfibo = [nfibx nfiby nfibz];
cur = 1; % disassemble into cell array again
for k = 1:length(tfibs),
    nfib{k} = nfibo(cur:(cur+size(tfibs{k},1)-1),:);
    cur = cur + size(tfibs{k},1);
end;


display('writing results');
t = inv(sdef.mat);
nfib = cellfun(@(x) x*t(1:3,1:3)' + repmat(t(1:3,4)',[size(x,1) 1]),nfib,'uniformoutput',false);

corrMy= diag(ones(1, 4)); corrMy(1:3, 4)= 1;  % switch to zero-based voxels
ftr.hMatrix =  q*sdef.mat *corrMy; 
ftr.curveSegCell = nfib;

return



function iy = load_def(name)
% Load a deformation field saved as an image

P      = [repmat(name,3,1), [',1,1';',1,2';',1,3']];
V      = spm_vol(P);
for k = 1:3,
    Def(:,:,:,k) = spm_load_float(V(k));
end;
mat    = V(1).mat;

iy.img = Def;

corrMy= diag(ones(1, 4)); corrMy(1:3, 4)= 1;  % switch to zero-based voxels
iy.mat = mat;% * corrMy;




return

