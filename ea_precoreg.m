function ea_precoreg(anchorImage, templateModality, preCoregImage, preCoregTransform)
% Pre-coregister pre-op anchor image to template (without reslicing).

% Make dirs
ea_mkdir(fileparts(preCoregImage));
ea_mkdir(fileparts(preCoregTransform));

% Copy anchor image to pre-coregistered image and leave anchor image untouched
if ~isfile(preCoregImage)
    copyfile(anchorImage, preCoregImage);
end

% Save transform
save_tmat = 1;

% Set flags for coregistration
flags.cost_fun = 'nmi';
flags.sep = [4 2];
flags.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
flags.fwhm = [7 7];

% Coregister and save transform
tmat = spm_coreg([ea_space,templateModality,'.nii,1'], [preCoregImage,',1'], flags);
tmat = spm_matrix(tmat(:)');
if save_tmat
    save(preCoregTransform, 'tmat');
end

% Set affine matrix to pre-coregistered image
preCoregVol = spm_vol(preCoregImage);
preCoregVol.mat = spm_get_space(preCoregImage, tmat \ preCoregVol.mat);
spm_write_vol(preCoregVol, spm_read_vols(preCoregVol));
