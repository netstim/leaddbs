function ea_create_dsistudio_src(dMRIfile,outputfilename)

[pth,fn,ext]=fileparts(outputfilename);
if strcmp(ext,'.gz')
   outputfilename=fullfile(pth,fn); 
end


nii=ea_load_nii(dMRIfile);

for grad=1:size(nii.img,4)
    tp=nii.img(:,:,:,grad);
    s.(['image',num2str(grad)])=tp(:);
end
s.voxel_size=nii.voxsize;
s.dimension=nii.dim;

[pth,fn,ext]=fileparts(dMRIfile);
bval=load(fullfile(pth,[fn,'.bval']));
bvec=load(fullfile(pth,[fn,'.bvec']));

s.b_table=[bval';bvec'];

save(outputfilename,'-struct','s','-v4');
gzip(outputfilename);
delete(outputfilename);
