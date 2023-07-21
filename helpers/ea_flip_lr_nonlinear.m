function XYZ=ea_flip_lr_nonlinear(from,to,interp)
% Flip files or coordinates (in mm) nonlinearly in the left-right direction
% based on pre-calculated asymmetric fliplr transformation.

ea_genflipspace; % will only perform if fliplr transformation doesn't exist

if ischar(from) % assume nifti file path
    if ~exist('interp','var')
        interp=4;
    end
    ea_flip_lr(from,to);
    options.subj.norm.log.method = fullfile(ea_space, 'fliplr', 'normmethod.json');
    options.subj.norm.transform.inverseBaseName = fullfile(ea_space, 'fliplr', 'InverseComposite.nii.gz');
    options.subj.norm.transform.forwardBaseName = fullfile(ea_space, 'fliplr', 'Composite.nii.gz');
    ea_apply_normalization_tofile(options,{to},{to},0,interp,to);
else % assume coordinate list
    spacedef = ea_getspacedef;
    ref = [ea_space, spacedef.templates{1}, '.nii'];
    from(:,1)=-from(:,1);

    % Map the 'mm' coords to 'vox' coords.
    XYZ_vox = ea_mm2vox(from, ref);
    
    % Flip the coordinates
    XYZ = ea_map_coords(XYZ_vox', ref, fullfile(ea_space, 'fliplr', 'Composite.nii.gz'), ref, 'ANTs');
    XYZ=XYZ';
end
