function varargout=ea_normalize_easyreg(options)

if ischar(options) % return name of method.
    varargout{1}='EasyReg (Iglesias 2023)';
    varargout{2}=1; % dummy output
    varargout{3}=0; % hassettings.
    varargout{4}=0; % is multispectral
    return
end

% Run EasyReg
[fwd_field, inv_field] = ea_easyreg([ea_space, options.primarytemplate, '.nii'], options.subj.coreg.anat.preop.(options.subj.AnchorModality));

% Move transformation file
ea_mkdir(fileparts(options.subj.norm.transform.forwardBaseName));
movefile(fwd_field, [options.subj.norm.transform.forwardBaseName, 'ants.nii.gz']);
movefile(inv_field, [options.subj.norm.transform.inverseBaseName, 'ants.nii.gz']);

% Apply registration
ea_apply_normalization(options)

%% add methods dump:
[scit, lcit] = ea_getspacedefcit;
cits = {'Iglesias, J. E. (2023). A ready-to-use machine learning tool for symmetric multi-modality registration of brain MRI. Scientific Reports, 13(1), Article 1. https://doi.org/10.1038/s41598-023-33781-0'};

if ~isempty(lcit)
    cits = [cits;{lcit}];
end

ea_methods(options,['Pre- (and post-) operative acquisitions were spatially normalized into ',ea_getspace,' space ',scit,' based on the preoperative acquisition (',options.subj.AnchorModality,') using the'...
    ' EasyReg approach as implemented in FreeSurfer (https://surfer.nmr.mgh.harvard.edu/fswiki/EasyReg).'],cits);
