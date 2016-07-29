function varargout=ea_normalize_spmnewseg(options)
% This is a function that normalizes both a copy of transversal and coronar
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
    switch spm('ver')
        case 'SPM12'
    varargout{1}='SPM12 Segment nonlinear [MR/CT]';
        case 'SPM8'
    varargout{1}='SPM8 New Segment nonlinear [MR/CT]';
    end
    varargout{2}={'SPM8','SPM12'};
    return
end

disp('This Normalization routine uses the advanced TPMs by Lorio 2016. See http://unil.ch/lren/home/menuinst/data--utilities.html');

if isfield(options.prefs, 'tranii_unnormalized')
    if exist([options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,'.gz'],'file')
        try
            gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,'.gz']);
        end
        try
            gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized,'.gz']);
        end
        try
            gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized,'.gz']);
        end
        try
            gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,'.gz']);
        end
    end
end

% First, do the coreg part:
try
    ea_coregmr(options,options.prefs.normalize.coreg);
end

% now segment the preoperative version.
disp('Segmenting preoperative version.');
ea_newseg([options.root,options.prefs.patientdir,filesep],options.prefs.prenii_unnormalized,0,options,0);
disp('*** Segmentation of preoperative MRI done.');

% Rename deformation fields:
movefile([options.root,options.patientname,filesep,'y_',options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,'y_ea_normparams.nii']);
movefile([options.root,options.patientname,filesep,'iy_',options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,'y_ea_inv_normparams.nii']);

% Apply estimated deformation to (coregistered) post-op data.

ea_apply_normalization(options)
