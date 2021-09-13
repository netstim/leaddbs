function template = ea_matchTemplate(moving, spacedef)
% Match moving image (using BIDS naming) to template to be used

if ~exist('spacedef', 'var')
    spacedef = ea_getspacedef;
end

parsedStruct = parseBIDSFilePath(moving);
movingModality = parsedStruct.suffix;

if ismember(movingModality, fieldnames(spacedef.norm_mapping))
    templateModality = spacedef.norm_mapping.(movingModality);
else
    templateModality = spacedef.misfit_template;
end

template = [ea_space, templateModality, '.nii'];
