function [itk_fwd_field, itk_inv_field] = ea_synthmorph(target_image, source_image)
    % Wrapper to run SynthMorph and get ANTs-like transforms

    %
    % SynthMorph
    %

    fs_fwd_field = strrep(source_image, '.nii', '_fs_fwd_field.nii');

    % Check Conda environment
    condaenv = ea_conda_env('SynthMorph');
    if ~condaenv.is_created
        ea_cprintf('CmdWinWarnings', 'Initializing SynthMorph conda environment...\n')
        condaenv.create;
        ea_cprintf('CmdWinWarnings', 'SynthMorph conda environment initialized.\n')
    elseif ~condaenv.is_up_to_date
        ea_cprintf('CmdWinWarnings', 'Updating SynthMorph conda environment...\n')
        condaenv.update;
        ea_cprintf('CmdWinWarnings', 'SynthMorph conda environment initialized.\n')
    end

    % Run SynthMorph
    synthmorph_exe = fullfile(ea_getearoot, 'ext_libs', 'SynthMorph', 'mri_synthmorph');
    synthmorph_cmd = {'python', ea_path_helper(synthmorph_exe), 'register', ...
        '-t', ea_path_helper(fs_fwd_field), ...
        ea_path_helper(source_image), ...
        ea_path_helper(target_image)};

    status = condaenv.system(strjoin(synthmorph_cmd, ' '));
    if status ~= 0
        ea_error('Registration using SynthMorph failed!', showdlg=false, simpleStack=true);
    end

    %
    % Convert transform
    %

    % Freesurfer to ITK transform (SynthMorph uses disp_ras format)
    itk_fwd_field = strrep(source_image, '.nii', '_itk_fwd_field.h5');
    ea_delete(itk_fwd_field);
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

    ea_delete(fs_fwd_field);

end


function [] = freesurfer_nii_to_itk_h5(warp_file_in, warp_file_out)

% substract mm coordinates for each voxel
n = load_nii(warp_file_in);
out = n.img;

% reshape output
out_rows = [-reshape(out(:,:,:,1),1,[]); -reshape(out(:,:,:,2),1,[]); reshape(out(:,:,:,3),1,[])];
out_column = reshape(out_rows,[],1);

% copy template h5 file
copyfile(fullfile(ea_getearoot, 'ext_libs', 'EasyReg', 'itk_h5_template.h5'), warp_file_out);

if ~strcmp(ea_getspace, 'MNI152NLin2009bAsym')
    % calculate TransformFixedParameters
    spacedef = ea_getspacedef;
    primarytemplate = [ea_space, spacedef.templates{1}, '.nii'];
    hdr = ea_fslhd(primarytemplate);
    TransformFixedParameters = zeros(18,1);
    TransformFixedParameters(1:3) = [hdr.dim1; hdr.dim2; hdr.dim3];
    TransformFixedParameters(4:6) = [-hdr.sto_xyz1(4); -hdr.sto_xyz2(4); hdr.sto_xyz3(4)]; % RAS to LPS applied
    TransformFixedParameters(7:9) = [hdr.pixdim1; hdr.pixdim2; hdr.pixdim3];
    TransformFixedParameters(10:18) = [-hdr.sto_xyz1(1:3)'/hdr.pixdim1; -hdr.sto_xyz2(1:3)'/hdr.pixdim2; hdr.sto_xyz3(1:3)'/hdr.pixdim3]; % RAS to LPS applied

    % update TransformFixedParameters in h5
    h5write(warp_file_out, "/TransformGroup/0/TransformFixedParameters", TransformFixedParameters);
end

% save TransformParameters in h5 
h5create(warp_file_out,"/TransformGroup/0/TransformParameters", numel(out_column));
h5write(warp_file_out, "/TransformGroup/0/TransformParameters", out_column);

end
