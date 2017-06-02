function XYZ=ea_flip_lr_nonlinear(from,to,interp)
% flip files / coords nonlinearly from Left to Right hemisphere based on asymmetric
% template
directory=[ea_space,'fliplr',filesep];

if ischar(from) % assume nifti file path
    if ~exist('interp','var')
        interp=4;
    end
    ea_genflipspace; % will only perform if doesnt exist
    options=ea_getptopts(directory);
    ea_apply_normalization_tofile(options,{from},{to},directory,0,interp,from);
    ea_flip_lr(to,to);
    tof=ea_load_nii(to);
    ea_reslice_nii(to,to,abs(tof.voxsize),[],[],[],[],[],0);
    
else % assume coordinate list
    spacedef=ea_getspacedef;
    % To map the 'mm' coords in src to the 'vox' coords in src:
    [~, XYZ_vox] = ea_map_coords([from,ones(size(from,1),1)]', [ea_space,spacedef.templates{1},'.nii']);
    XYZ = ea_map_coords(XYZ_vox, [ea_space,spacedef.templates{1},'.nii'], [directory,'y_ea_inv_normparams.nii'], '');
    XYZ(1)=-XYZ(1);
    XYZ=XYZ';
end