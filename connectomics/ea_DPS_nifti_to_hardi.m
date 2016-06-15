function [mrs nii] = ea_DPS_nifti_to_hardi(niiName,bvecname,bvalname)

nii_orig = ea_load_nii(niiName);
nii=nii_orig(1);
nii.img=zeros([size(nii.img),length(nii_orig)]);

for grad=1:length(nii_orig);
   nii.img(:,:,:,grad)=nii_orig(grad).img; 
end


mrs = ea_DPS_nifti_to_mrs(nii);

if nargin == 1,
    bvec = load([niiName(1:end-3), 'bvec']);
    bval = load([niiName(1:end-3), 'bval']);
else
    bvec = load(bvecname);
    bval = load(bvalname);    
end

if size(bvec,1)>size(bvec,2)
    bvec=bvec';
end
if size(bval,1)>size(bval,2)
    bval=bval';
end
T = 1;

for k = 1:size(bval,2),
    gdir = T*bvec(:,k);
    gdir = gdir / (eps+norm(gdir));
    tensor(:,:,k) = gdir*gdir' *bval(k);
end;

mrs.user.bfactor = bval;
mrs.user.bDir = bvec;
mrs.user.bTensor = tensor;


