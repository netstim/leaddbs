function varargout=ea_normalize_spmnewseg(options)
% This is a function that normalizes both a copy of transversal and coronal
% images into MNI-space. The goal was to make the procedure both robust and
% automatic, but still, it must be said that normalization results should
% be taken with much care because all reconstruction results heavily depend
% on these results. Normalization of DBS-MR-images is especially
% problematic since usually, the field of view doesn't cover the whole
% brain (to reduce SAR-levels during acquisition) and since electrode
% artifacts can impair the normalization process. Therefore, normalization
% might be best archieved with other tools that have specialized on
% normalization of such image data.
%
% The procedure used here uses the SPM8 "New Segment", or SPM12 "Segment" routine and
% is probably the most straight-forward way using SPM8.
%
% This function uses resize_img.m authored by Ged Rigdway
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='SPM12 Segment (Ashburner 2005)';
    varargout{2}=1; % dummy output
    varargout{3}=1; % hassettings.
    varargout{4}=1; % is multispectral
    return
end

fprintf('\nThis Normalization routine uses the advanced TPMs by Lorio 2016.\nSee http://unil.ch/lren/home/menuinst/data--utilities.html\n\n');

% Run SPM New Segment
disp('Segmenting preoperative version...');
preopImages = struct2cell(options.subj.coreg.anat.preop);
ea_newseg(preopImages, 0, 0, 1);
disp('Segmentation of preoperative MRI done.');

[directory, preopAnchorName] = fileparts(preopImages{1});
directory = [directory, filesep];
preopAnchorName = [preopAnchorName, '.nii'];

disp('done.');

% Rename Segmentations (c1, c2, c3)
mod = replace(options.subj.AnchorModality, textBoundary('start') + alphanumericsPattern + "_", "");
movefile([directory, 'c1', preopAnchorName], setBIDSEntity(preopImages{1}, 'mod', mod, 'label', 'GM', 'suffix', 'mask'));
movefile([directory, 'c2', preopAnchorName], setBIDSEntity(preopImages{1}, 'mod', mod, 'label', 'WM', 'suffix', 'mask'));
movefile([directory, 'c3', preopAnchorName], setBIDSEntity(preopImages{1}, 'mod', mod, 'label', 'CSF', 'suffix', 'mask'));

% Rename deformation fields
ea_mkdir(fileparts(options.subj.norm.transform.forwardBaseName));
movefile([directory, 'y_', preopAnchorName], [options.subj.norm.transform.forwardBaseName, 'spm.nii']);
movefile([directory, 'iy_', preopAnchorName], [options.subj.norm.transform.inverseBaseName, 'spm.nii']);

% Apply estimated deformation to (coregistered) post-op images.
ea_apply_normalization(options)

modality = regexp(preopImages, '(?<=_)[^\W_]+(?=\.nii(\.gz)?$)', 'match', 'once');

ea_methods(options,['Pre- (and post-) operative acquisitions were spatially normalized into ',ea_getspace,' space based on preoperative acquisition(s) (',strjoin(modality, ', '),') using the'...
    'Unified Segmentation Approach as implemented in SPM12 (Ashburner 2005; www.fil.ion.ucl.ac.uk/spm/software/spm12/).'], ...
    {'Ashburner, J., & Friston, K. J. (2005). Unified segmentation., 26(3), 839?851. http://doi.org/10.1016/j.neuroimage.2005.02.018'});
