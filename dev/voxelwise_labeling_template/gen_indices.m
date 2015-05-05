% this voxel-wise parcellation scheme has been built using the
% mni_icbm152_gm_tal_nlin_asym_09c template. The template was downsampled
% to a voxel-size of 3.2 mm isotropic and thresholded at 0.4 (threshold was
% actually set to 0.40112 to obtain a round number of exactly 35k nodes).
% Indices start on the left hemisphere from dorsal to ventral and then head
% over to the right hemisphere (also dorsal to ventral).

clear
clc
V=spm_vol('mni_gm_mask_t04.nii');
X=spm_read_vols(V);
V.fname='mni_parcellation_t04.nii';

[xx,yy,zz]=ind2sub(size(X),find(X(:)));
XYZ_vx=[xx,yy,zz,ones(length(xx),1)]';
XYZ_mm=V.mat*XYZ_vx;
XYZ_mm=XYZ_mm(1:3,:);
XYZ_vx=XYZ_vx(1:3,:);
xxes=XYZ_mm(1,:);
[posix]=find(xxes>0);
[negix]=find(xxes<0);

sequence=zeros(length(xxes),1);
cnt=1;
for entry=1:length(posix);
    sequence(cnt)=posix(entry);
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
of=fopen('mni_parcellation_t04.txt','w');
for entry=1:length(XYZ_vx)
   X(XYZ_vx(1,entry),XYZ_vx(2,entry),XYZ_vx(3,entry))=cnt;
fprintf(of,'%d %s \n',entry,[num2str(XYZ_mm(1,entry)),'_',num2str(XYZ_mm(2,entry)),'_',num2str(XYZ_mm(3,entry))]);
   cnt=cnt+1;
end
spm_write_vol(V,X);
fclose(of);
of=fopen('mni_parcellation_t04.txt');

atlas_lgnd=textscan(of,'%d %s');

