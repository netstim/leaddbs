% this voxel-wise parcellation scheme has been built using the
% mni_icbm152_gm_tal_nlin_asym_09c template. The template was downsampled
% to a voxel-size of 5 mm isotropic and thresholded at 0.5 (threshold was
% actually set to 0.5 to obtain a round number of 8057 nodes).
% Indices start on the left hemisphere from dorsal to ventral and then head
% over to the right hemisphere (also dorsal to ventral).

clear
clc
ea_reslice_nii([ea_getearoot,'templates/mni_icbm152_gm_tal_nlin_asym_09c.nii'],'mni_gm_05.nii',[5,5,5]);
V=spm_vol('mni_gm_05.nii');
X=spm_read_vols(V);
X(X<0.5)=0;
X=logical(X);
disp(['Sum is ',num2str(sum(X(:))),'.']);

V.fname='tmni_gm_05.nii';
V.dt(1) = 16;
spm_write_vol(V,X);


V=spm_vol('tmni_gm_05.nii');
X=spm_read_vols(V);
disp(['Sum is ',num2str(sum(X(:))),'.']);

V.fname='mni_parcellation_05_t05.nii';

[xx,yy,zz]=ind2sub(size(X),find(X(:)));
XYZ_vx=[xx,yy,zz,ones(length(xx),1)]';
XYZ_mm=V.mat*XYZ_vx;
XYZ_mm=XYZ_mm(1:3,:);
XYZ_vx=XYZ_vx(1:3,:);
xxes=XYZ_mm(1,:);
[posix]=find(xxes>0);
[negix]=find(xxes<0);
zeroix=find(xxes==0);

sequence=zeros(length(xxes),1);
cnt=1;
for entry=1:length(posix);
    sequence(cnt)=posix(entry);
    cnt=cnt+1;
end
for entry=1:length(zeroix);
    sequence(cnt)=zeroix(entry);
    cnt=cnt+1;
end
for entry=1:length(negix);
    sequence(cnt)=negix(entry);
    cnt=cnt+1;
end

sequence=flipud(sequence);

XYZ_vx=XYZ_vx(:,sequence);
XYZ_mm=XYZ_mm(:,sequence);
cnt=1;
X(:)=0;
of=fopen('mni_parcellation_05_t05.txt','w');
for entry=1:length(XYZ_vx)
   X(XYZ_vx(1,entry),XYZ_vx(2,entry),XYZ_vx(3,entry))=cnt;
fprintf(of,'%d %s \n',entry,[num2str(XYZ_mm(1,entry)),'_',num2str(XYZ_mm(2,entry)),'_',num2str(XYZ_mm(3,entry))]);
   cnt=cnt+1;
end
V.dt(1) = 8;
spm_write_vol(V,X);
fclose(of);
of=fopen('mni_parcellation_05_t05.txt');

atlas_lgnd=textscan(of,'%d %s');

