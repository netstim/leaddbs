function XYZ=ea_flip_lr_nonlinear(from,to,interp)
% flip files / coords nonlinearly from Left to Right hemisphere based on asymmetric
% template

directory=[ea_space,'fliplr',filesep];
ea_genflipspace; % will only perform if doesnt exist
if ischar(from) % assume nifti file path
    if ~exist('interp','var')
        interp=4;
    end
    options=ea_getptopts(directory);
    ea_flip_lr(from,to);
    ea_apply_normalization_tofile(options,{to},{to},0,interp,to);
else % assume coordinate list
    spacedef=ea_getspacedef;
    from(:,1)=-from(:,1);

    % To map the 'mm' coords in src to the 'vox' coords in src:
    [~, XYZ_vox] = ea_map_coords([from,ones(size(from,1),1)]', [ea_space,spacedef.templates{1},'.nii']);
    XYZ = ea_map_coords(XYZ_vox, [ea_space,spacedef.templates{1},'.nii'], [directory,'inverseTransform'], '');
    XYZ=XYZ';
end
