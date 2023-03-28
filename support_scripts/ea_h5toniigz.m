function ea_h5toniigz(h5file,output_file,reference_file)


if ~exist('reference_file','var') || isempty(reference_file)
    if contains(output_file,'to-anchorNative')
        anchorNative_files = dir(fullfile(fileparts(fileparts(fileparts(output_file))), 'coregistration', 'anat', '*space-anchorNative*.nii'));
        if isempty(anchorNative_files)
            error('No coregistered images found to use as reference')
        end
        reference_file = fullfile(anchorNative_files(1).folder, anchorNative_files(1).name)
    else
        reference_file = fullfile(ea_space,'t1.nii');
    end
end

ea_libs_helper;
basedir = [fullfile(ea_getearoot,'ext_libs','ANTs'), filesep];
if ispc
    applyTransforms = ea_path_helper([basedir, 'antsApplyTransforms.exe']);
else
    applyTransforms = ea_path_helper([basedir, 'antsApplyTransforms.', computer('arch')]);
end

%antsApplyTransforms -t glanatInverseComposite.h5 -r anat_t2.nii -o [glanatInverseComposite.nii.gz,1] --float -v 1
cmd = [applyTransforms, ...
       ' -t ', h5file ...
       ' -r ', reference_file ...
       ' -o [', output_file ...
       ',1]' ...
       ' --float 1'...
       ' -v 1'];
if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end