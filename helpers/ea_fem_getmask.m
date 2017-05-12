function fv=ea_fem_getmask(options)

% load nifti
try
    switch options.native
        case 0 % template space
            nii=ea_load_nii([ea_space(options,'atlases'),options.atlasset,filesep,'gm_mask.nii']);
        case 1 % native space
            nii=ea_load_nii([options.root,options.patientname,filesep,'atlases',filesep,options.atlasset,filesep,'gm_mask.nii']);
    end
catch
    ea_error('The selected atlas set seems incompatible with this approach.');
end

[xx,yy,zz]=ind2sub(size(nii.img),find(nii.img>0)); %(mean(nii.img(nii.img~=0))/3))); % find 3D-points that have correct value.
if isempty(xx)
    ea_error('The selected atlas set seems incompatible with this approach.');
end
XYZ=[xx,yy,zz]; % concatenate points to one matrix.
XYZ=nii.mat*[XYZ,ones(length(xx),1)]'; % map to mm-space
XYZ=XYZ(1:3,:)';

bb=[0,0,0,1;size(nii.img),1];

bb=nii.mat*bb';
bb=bb(1:3,:)';
gv=cell(3,1);
for dim=1:3
    gv{dim}=linspace(bb(1,dim),bb(2,dim),size(nii.img,dim));
end
[X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});
if options.prefs.hullsmooth
    nii.img = smooth3(nii.img,'gaussian',options.prefs.hullsmooth);
end
fv=isosurface(X,Y,Z,permute(nii.img,[2,1,3]),max(nii.img(:))/2);
fvc=isocaps(X,Y,Z,permute(nii.img,[2,1,3]),max(nii.img(:))/2);
fv.faces=[fv.faces;fvc.faces+size(fv.vertices,1)];
fv.vertices=[fv.vertices;fvc.vertices];




%figure
%patch('vertices',fv.vertices,'faces',fv.faces,'facecolor','r');
