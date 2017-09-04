function ea_gendistheatmap(M,regno,pts,fovimg)
if ~exist('fovimg','var')
    fovimg=[ea_space,'atlas.nii'];
end
if ~exist('regno','var') || isempty(regno)
    regno=1;
end
if ~exist('pts','var') || isempty(pts)
    pts=1:length(M.patient.list);
end

fovimg=ea_load_nii(fovimg);

fovimg.fname=[M.root,M.clinical.labels{regno},'_dist_heatmap.nii'];

[xx,yy,zz]=ind2sub(size(fovimg.img),1:numel(fovimg.img));
XYZ=[xx;yy;zz;ones(1,length(xx))];
XYZ=fovimg.mat*XYZ;
XYZ=XYZ(1:3,:)';

for pt=1:length(M.patient.list)
    acs(pt,:)=mean([mean(M.elstruct(pt).coords_mm{1}(logical(M.S(pt).activecontacts{1}(1:4)),:),1);...
        mean(M.elstruct(pt).coords_mm{2}(logical(M.S(pt).activecontacts{2}(1:4)),:),1).*[-1,1,1]]);
end

% figure, plot3(acs(:,1),acs(:,2),acs(:,3),'r*');
% axis equal
ea_dispercent(0,'Iterating voxels');
dimen=length(XYZ);
N=length(M.patient.list(pts));
I=M.clinical.vars{regno}(pts);
for vx=1:dimen
   %D=squareform(pdist([XYZ(vx,:);acs]));
   D=-pdist([XYZ(vx,:);acs]);
   D=D(1:N)';
   fovimg.img(vx)=corr(D,I);
   ea_dispercent(vx/dimen);
end
ea_dispercent(1,'end');


fovimg.dt=[16,0];
ea_write_nii(fovimg);

