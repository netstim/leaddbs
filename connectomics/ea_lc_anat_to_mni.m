function [mm, vox] = ea_lc_anat_to_mni(options, points, anatomy, template)
% Map from the recorded B0 target frame into the active normalization frame.
anchor = options.subj.coreg.anat.preop.(options.subj.AnchorModality);
preproc = options.subj.preopAnat.(options.subj.AnchorModality).preproc;
if strcmp(anatomy, anchor)
    anchorPoints = points;
elseif strcmp(anatomy, preproc)
    % ea_precoreg changes the anchor's world frame before normalization.
    original = ea_get_affine(anatomy);
    target = ea_get_affine(anchor);
    if norm(original-target, 'fro') < 1e-4
        bridge = eye(4);
    else
        file = options.subj.coreg.transform.(options.subj.AnchorModality);
        if ~ischar(file) || ~isfile(file)
            error('LeadDBS:MissingPrecoreg', 'The recorded anatomical frame requires its pre-coregistration transform.');
        end
        saved = load(file, 'tmat');
        bridge = inv(saved.tmat);
        if norm(bridge*original-target, 'fro') > 1e-3
            error('LeadDBS:InconsistentPrecoreg', 'The saved pre-coregistration transform does not match the anatomical headers.');
        end
    end
    mapped = target \ (bridge * original * [points; ones(1,size(points,2))]);
    anchorPoints = mapped(1:3,:);
else
    error('LeadDBS:UnknownFiberAnatomy', 'Fiber anatomy must be the protocol anchor or its preprocessing image.');
end
% Use the standard subject-aware path: SPM is converted to ITK as needed,
% ANTs point direction is handled centrally, and FSL keeps its own mapping.
[mm, vox] = ea_map_coords(anchorPoints, anchor, ...
    fullfile(options.root, options.patientname, 'inverseTransform'), template);
end

