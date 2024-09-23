function varargout=ea_normalize_synthmorph(options)

if ischar(options) % return name of method.
    varargout{1}='SynthMorph (Hoffmann 2024)';
    varargout{2}=1; % dummy output
    varargout{3}=0; % hassettings.
    varargout{4}=0; % is multispectral
    return
end

% Run SynthMorph
[fwd_field, inv_field] = ea_synthmorph([ea_space, options.primarytemplate, '.nii'], options.subj.coreg.anat.preop.(options.subj.AnchorModality));

% Move transformation file
ea_mkdir(fileparts(options.subj.norm.transform.forwardBaseName));
movefile(fwd_field, [options.subj.norm.transform.forwardBaseName, 'ants.nii.gz']);
movefile(inv_field, [options.subj.norm.transform.inverseBaseName, 'ants.nii.gz']);

% Apply registration
ea_apply_normalization(options)

%% add methods dump:
[scit, lcit] = ea_getspacedefcit;
cits = {'Hoffmann, M., Billot, B., Greve, D.N., Iglesias, J.E., Fischl, B., Dalca, A.V., 2022. SynthMorph: Learning contrast-invariant registration without acquired images. IEEE Trans. Med. Imag. 41, 543–558. https://doi.org/10.1109/TMI.2021.3116879'
        'Hoffmann, M., Hoopes, A., Fischl, B., Dalca, A.V., 2023. Anatomy-specific acquisition-agnostic affine registration learned from fictitious images, in: Medical Imaging 2023: Image Processing. Presented at the Medical Imaging 2023: Image Processing, SPIE, p. 1246402. https://doi.org/10.1117/12.2653251'
        'Hoffmann, M., Hoopes, A., Greve, D.N., Fischl, B., Dalca, A.V., 2024. Anatomy-aware and acquisition-agnostic joint registration with SynthMorph. Imaging Neuroscience 2, 1–33. https://doi.org/10.1162/imag_a_00197'};
if ~isempty(lcit)
    cits = [cits;{lcit}];
end

ea_methods(options,['Pre- (and post-) operative acquisitions were spatially normalized into ',ea_getspace,' space ',scit,' based on the preoperative acquisition (',options.subj.AnchorModality,') using the'...
    ' SynthMorph approach as implemented in FreeSurfer (https://synthmorph.io).'],cits);
