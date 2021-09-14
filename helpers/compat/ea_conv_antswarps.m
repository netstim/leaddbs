function ea_conv_antswarps(transform_file_name, reference, float)
% switches ants transform extenstion between .h5 and .nii.gz

[~,~,ext] = fileparts(transform_file_name);

switch ext
    case '.gz'
        ext = '.nii.gz';
        out_ext = '.h5';
    case '.h5'
        out_ext = '.nii.gz';
    otherwise
        error(['Unrecognized ANTs transform file extension: ' ext])
end

out_file_name = strrep(transform_file_name, ext, out_ext);

antsdir=[ea_getearoot,'ext_libs',filesep,'ANTs',filesep];
if ispc
    applyTransforms = ea_path_helper([antsdir, 'antsApplyTransforms.exe']);
else
    applyTransforms = [antsdir, 'antsApplyTransforms.', computer('arch')];
end

cmd = [applyTransforms ' -r ' reference ' -t ' transform_file_name ' -o [' out_file_name ',1]'];

if exist('float', 'var')
    if ischar(float) && strcmp(float, 'float') || float
        cmd = [cmd, ' --float'];
    end
end

cmd = [cmd, ' -v 1'];

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end

ea_delete(transform_file_name)

