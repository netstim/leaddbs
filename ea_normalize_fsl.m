function varargout=ea_normalize_fsl(options)
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
% The procedure used here uses the ANTs Syn approach to map a patient's
% brain to MNI space directly.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


if ischar(options) % return name of method.
    varargout{1}='FNIRT (Andersson 2010)';
    varargout{2}=1; % dummy output
    varargout{3}=1; % hassettings.
    varargout{4}=0; % is multispectral
    return
end

% FSL FNIRT nonlinear registration
ea_fnirt([ea_space, options.primarytemplate, '.nii'],...
         options.subj.coreg.anat.preop.(options.subj.AnchorModality),...
         options.subj.norm.anat.preop.(options.subj.AnchorModality));

% Clean up FNIRT log file
logFileBase = ea_niifileparts(options.subj.coreg.anat.preop.(options.subj.AnchorModality));
ea_delete([logFileBase, '_to_*.log']);

% Move transformation file
transformBase = ea_niifileparts(options.subj.norm.anat.preop.(options.subj.AnchorModality));
transform = dir([transformBase, 'WarpField*']);
ext = regexp(transform(end).name, '(?<=WarpField)\.nii(\.gz)?$', 'match', 'once');
movefile([transformBase, 'WarpField', ext], [options.subj.norm.transform.forwardBaseName, 'fnirt', ext]);
movefile([transformBase, 'InverseWarpField', ext], [options.subj.norm.transform.inverseBaseName, 'fnirt', ext]);

% Apply registration
ea_apply_normalization(options)

%% add methods dump:
[scit, lcit] = ea_getspacedefcit;
cits = {
    'Andersson JLR, Jenkinson M, Smith S (2010) Non-linear registration, aka spatial normalisation. FMRIB technical report TR07JA2'
    'Woolrich MW, Jbabdi S, Patenaude B, Chappell M, Makni S, Behrens T, Beckmann C, Jenkinson M, Smith SM (2009) Bayesian analysis of neuroimaging data in FSL. NeuroImage, 45:S173-86'
    };

if ~isempty(lcit)
    cits = [cits;{lcit}];
end

ea_methods(options,['Pre- (and post-) operative acquisitions were spatially normalized into ',ea_getspace,' space ',scit,' based on the preoperative acquisition (',options.subj.AnchorModality,') using the'...
    ' FNIRT approach as implemented in the FMRIB Software Library (Andersson 2010; Woolrich 2009 https://fsl.fmrib.ox.ac.uk/).'],cits);
