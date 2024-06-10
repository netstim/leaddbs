function [itk_fwd_field, itk_inv_field] = ea_easyreg(target_image, source_image)
    % Wrapper to run EasyReg and get ANTs-like transforms

    %
    % EasyReg
    %

    target_seg = [target_image(1:end-4) '_synthseg.nii'];
    source_seg = [source_image(1:end-4) '_synthseg.nii'];
    fs_fwd_field = [source_image(1:end-4) '_fs_fwd_field.nii'];

    % Check Conda environment
    condaenv = ea_conda_env('EasyReg');
    if ~condaenv.is_created
        ea_cprintf('CmdWinWarnings', 'Initializing EasyReg conda environment...\n')
        condaenv.create;
        ea_cprintf('CmdWinWarnings', 'EasyReg conda environment initialized.\n')
    elseif ~condaenv.is_up_to_date
        ea_cprintf('CmdWinWarnings', 'Updating EasyReg conda environment...\n')
        condaenv.update;
        ea_cprintf('CmdWinWarnings', 'EasyReg conda environment initialized.\n')
    end

    % Run EasyReg
    easyreg_exe = fullfile(ea_getearoot, 'ext_libs', 'EasyReg', 'mri_easyreg');
    easyreg_cmd = {'python', ea_path_helper(easyreg_exe), ...
        '--ref', ea_path_helper(target_image), '--ref_seg', ea_path_helper(target_seg), ...
        '--flo', ea_path_helper(source_image), '--flo_seg', ea_path_helper(source_seg), ...
        '--fwd_field', ea_path_helper(fs_fwd_field), ...
        '--threads -1'};
    status = condaenv.system(strjoin(easyreg_cmd, ' '));
    if status ~= 0
        ea_error('Registration using EasyReg failed!', showdlg=false, simpleStack=true);
    end

    %
    % Convert transform
    %

    % Freesurfer to ITK transform
    itk_fwd_field = [source_image(1:end-4) '_itk_fwd_field.h5'];
    freesurfer_nii_to_itk_h5(fs_fwd_field, itk_fwd_field);

    % Set-up Custom Slicer
    s4l = ea_slicer_for_lead;
    if ~s4l.is_up_to_date()
        s4l.install();
    end

    % Invert transform
    itk_inv_field = strrep(itk_fwd_field, '_fwd_', '_inv_');
    ea_slicer_invert_transform(itk_fwd_field, source_image, itk_inv_field)

    % .h5 to .nii.gz
    ea_conv_antswarps(itk_fwd_field, target_image, 1);
    ea_conv_antswarps(itk_inv_field, source_image, 1);

    itk_fwd_field = strrep(itk_fwd_field, '.h5', '.nii.gz');
    itk_inv_field = strrep(itk_inv_field, '.h5', '.nii.gz');

    ea_delete({source_seg, fs_fwd_field});

end


function [] = freesurfer_nii_to_itk_h5(warp_file_in, warp_file_out)

% substract mm coordinates for each voxel

n = load_nii(warp_file_in);
s = n.hdr.dime.dim(2:4);
index = 1:prod(s);
[v1,v2,v3] = ind2sub(s,index);
mm = ea_vox2mm([v1',v2',v3'], get_mat);
mm = reshape(mm, [s,3]);
out = n.img - mm;

% reshape output

out_rows = [-reshape(out(:,:,:,1),1,[]); -reshape(out(:,:,:,2),1,[]); reshape(out(:,:,:,3),1,[])];
out_column = reshape(out_rows,[],1);

% save

copyfile(fullfile(ea_getearoot, 'ext_libs', 'EasyReg', 'itk_h5_template.h5'), warp_file_out);
h5create(warp_file_out,"/TransformGroup/0/TransformParameters", numel(out_column));
h5write(warp_file_out,"/TransformGroup/0/TransformParameters", out_column);

end

function mat = get_mat()

mat = [ 0.5	0	0	 -98.5;...
        0	0.5	0	-134.5;...
        0	0	0.5	 -72.5;...
        0	0	0	   1];

end