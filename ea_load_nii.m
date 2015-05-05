function nii=ea_load_nii(fname)
% simple nifti reader using SPM.

if strcmp(fname(end-2:end),'.gz')
    wasgzip=1;
    gunzip(fname);
    fname=fname(1:end-3);
else
    wasgzip=0;
end

    nii=spm_vol(fname);
    
img=spm_read_vols(nii);
if length(nii)>1 % multi volume;
   for n=1:length(nii)
      nii(n).img=squeeze(img(:,:,:,n));
   end
else
nii.img=img;
end


try
nii.hdr.dime.pixdim=nii.mat(logical(eye(4)));
end
if wasgzip
    delete(fname); % since gunzip makes a copy of the zipped file.
end