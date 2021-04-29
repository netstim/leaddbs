function fv=ea_nii2fv_v2(nii_filename,smooth,smooth_value)

nii = ea_load_nii(nii_filename);
nii.img(isnan(nii.img))=0;

bb=[1,1,1;size(nii.img)];
bb=ea_vox2mm(bb,nii.mat);
gv=cell(3,1);
for dim=1:3
    gv{dim}=linspace(bb(1,dim),bb(2,dim),size(nii.img,dim));
end
[X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});

if smooth
    nii.img = smooth3(nii.img,'gaussian',smooth_value);
    %options.prefs.hullsmooth);
end
%toc

fv=isosurface(X,Y,Z,permute(nii.img,[2,1,3]),max(nii.img(:))/2);
fvc=isocaps(X,Y,Z,permute(nii.img,[2,1,3]),max(nii.img(:))/2);
fv.faces=[fv.faces;fvc.faces+size(fv.vertices,1)];
fv.vertices=[fv.vertices;fvc.vertices];
