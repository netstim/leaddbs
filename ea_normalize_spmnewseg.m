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

directory=[options.root,options.patientname,filesep];
if isfield(options.prefs, 'tranii_unnormalized')
    if exist([directory,options.prefs.tranii_unnormalized,'.gz'],'file')
        try
            gunzip([directory,options.prefs.tranii_unnormalized,'.gz']);
        end
        try
            gunzip([directory,options.prefs.cornii_unnormalized,'.gz']);
        end
        try
            gunzip([directory,options.prefs.sagnii_unnormalized,'.gz']);
        end
        try
            gunzip([directory,options.prefs.prenii_unnormalized,'.gz']);
        end
    end
end

% now segment the preoperative version.
disp('Segmenting preoperative version...');

ea_newseg_pt(options,0,0,1);
disp('done.');

% Rename deformation fields:
try copyfile([directory,'y_',options.prefs.prenii_unnormalized],[directory,'y_ea_normparams.nii']); end
try copyfile([directory,'iy_',options.prefs.prenii_unnormalized],[directory,'y_ea_inv_normparams.nii']); end

% Apply estimated deformation to (coregistered) post-op data.

ea_apply_normalization(options)

[~,anatpresent]=ea_assignpretra(options);

ea_methods(options,['Pre- (and post-) operative acquisitions were spatially normalized into ',ea_getspace,' space based on preoperative acquisition(s) (',ea_cell2strlist(anatpresent),') using the'...
    'Unified Segmentation Approach as implemented in SPM12 (Ashburner 2005; www.fil.ion.ucl.ac.uk/spm/software/spm12/).',...
    ],...
    {...
    'Ashburner, J., & Friston, K. J. (2005). Unified segmentation., 26(3), 839?851. http://doi.org/10.1016/j.neuroimage.2005.02.018',...
    });
