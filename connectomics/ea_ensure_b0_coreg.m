function options = ea_ensure_b0_coreg(options)
% Use the current registration contract; do not discover transforms.
directory = fullfile(options.root, options.patientname);
protocol = fullfile(directory, 'coregistration', 'log', 'ea_lc_b0_coreg.mat');
if isfile(protocol)
    ea_lc_coreg_record(options);
    return;
end
[options, ~, b0] = ea_lc_resolve_b0(options);
[options, anatRel, anatName] = ea_lc_resolve_anat_anchor(options);
if isempty(b0) || isempty(anatRel)
    error('LeadDBS:MissingB0CoregInput', 'B0 and anatomical anchor are required for coregistration.');
end
[~, b0Name] = ea_niifileparts(b0);
output = fullfile(directory, 'preprocessing', 'anat', [b0Name, '2', anatName, '.nii']);
ea_coregimages(options, b0, fullfile(directory, anatRel), output, {}, 1, [], 1);
ea_lc_coreg_record(options);
end
