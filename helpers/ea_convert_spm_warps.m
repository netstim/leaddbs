function ea_convert_spm_warps(subj)
% Convert SPM deformation field into ITK format

forwardBaseName = subj.norm.transform.forwardBaseName;
inverseBaseName = subj.norm.transform.inverseBaseName;

if isfile([forwardBaseName, 'spm.nii'])
    ea_spm_fwd_displacement_field_to_ants([forwardBaseName, 'spm.nii'], [forwardBaseName, 'ants.nii.gz']);
    ea_slicer_invert_transform([forwardBaseName, 'ants.nii.gz'], subj.coreg.anat.preop.(subj.AnchorModality), [inverseBaseName, 'ants.nii.gz']);

    ea_delete({[forwardBaseName, 'spm.nii']; [inverseBaseName, 'spm.nii']});
end
