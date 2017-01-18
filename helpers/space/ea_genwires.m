function ea_genwires
load([ea_space,'ea_space_def.mat'])

nii=ea_load_nii([ea_space,spacedef.templates{1},'.nii']);
wires=ea_detect_edges_3d(nii.img,0.2);
wires=wires>0.1;
wires=wires*255;
%figure, imagesc(squeeze(wires(round(size(wires,1)/2),:,:)));
save([ea_space,'wires.mat'],'wires');
