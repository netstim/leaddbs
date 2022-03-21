function ea_resliceanat(preopAnchorImage)
% Reslice the pre-op anchor image if the resolution is not sufficient

V = spm_vol(preopAnchorImage);
dim = V.mat(logical(eye(4)));
dim = abs(dim(1:3));

if any(dim>0.7)
    fprintf('\nInterpolating pre-op anchor image to increase resolution to 0.7...')
    ea_resample_image_by_spacing(preopAnchorImage,[0.7 0.7 0.7],0,5,0,preopAnchorImage);
    % ea_reslice_nii(preopAnchorImage,preopAnchorImage,[0.7 0.7 0.7],0,[],1);
    fprintf('\n\n');
end
