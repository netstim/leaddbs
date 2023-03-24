function ea_h5toniigz(h5file,output_file)
ea_libs_helper;
basedir = [fullfile(ea_getearoot,'ext_libs','ANTs'), filesep];
if ispc
    applyTransforms = ea_path_helper([basedir, 'antsApplyTransforms.exe']);
else
    applyTransforms = ea_path_helper([basedir, 'antsApplyTransforms.', computer('arch')]);
end
mnifile = fullfile(ea_space,'t1.nii');
%antsApplyTransforms -t glanatInverseComposite.h5 -r anat_t2.nii -o [glanatInverseComposite.nii.gz,1] --float -v 1
cmd = [applyTransforms, ...
       ' -t ', h5file ...
       ' -r ', mnifile ...
       ' -o [', output_file ...
       ',1]' ...
       ' --float 1'...
       ' -v 1'];
if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end