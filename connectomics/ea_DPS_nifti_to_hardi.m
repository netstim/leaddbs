function [mrs nii_orig] = DPS_nifti_to_hardi(niiName,bvecname,bvalname)

nii_orig = load_untouch_nii(niiName);
nii = nii_orig;

mrs = DPS_nifti_to_mrs(nii);

if nargin == 1,
    bvec = importdata([niiName(1:end-3), 'bvec']);
    bval = importdata([niiName(1:end-3), 'bval']);
else
    bvec = importdata(bvecname);
    bval = importdata(bvalname);    
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

end

