function ea_get_field_from_csv(ref_image, Field_array_file, Activation_threshold_VTA, sideLabel, outputBasePath)

ea_dispt('Creating nifti header for export...');
% create nifti
[~, ~, endian] = computer;
switch endian
    case 'L'
        endian = 0;
    case 'B'
        endian = 1;
end

% split to corrdinates and field
Field_array = table2array(readtable(Field_array_file));
Field_coords = Field_array(:,2:4);
Field_vals = Field_array(:,8);  % others are the components

% convert to native voxel space (will be as floating numbers)
Field_vox_native = ea_mm2vox(Field_coords, ref_image)';
ref_image_nii = ea_load_nii(ref_image);
array_VTA = zeros(size(ref_image_nii.img));

% there will be overwriting due to lower resolution
for point_i = 1:size(Field_vox_native,2)
    x_ind = round(Field_vox_native(1,point_i));
    y_ind = round(Field_vox_native(2,point_i));
    z_ind = round(Field_vox_native(3,point_i));
    array_VTA(x_ind, y_ind,z_ind) = Field_vals(point_i);
    %disp(array_VTA(x_ind, y_ind,z_ind))
end

Field_ref.mat = ref_image_nii.mat;
Field_ref.dim=size(ref_image_nii.img);
Field_ref.dt = [4, endian];
Field_ref.n=[1 1];
Field_ref.descrip='oss-dbs-v2 - Field_ref';

[filepath,filename,~] = fileparts(Field_array_file);

Field_ref.img = array_VTA; 
Field_ref.img(isnan(Field_ref.img)) = 0;
Field_ref.img(Field_ref.img>1000.0) = 1000.0;
Field_ref.fname = [outputBasePath, 'efield_model-ossdbs_hemi-', sideLabel, '.nii'];
ea_write_nii(Field_ref);
%ea_autocrop([outputBasePath, 'efield_model-ossdbs_hemi-', sideLabel, '.nii']);

% also create VATs directly
VTA_interp = array_VTA >= (Activation_threshold_VTA);
Vvat2 = Field_ref;
Vvat2.descrip='oss-dbs-v2 - VAT_ref';
Vvat2.fname = [outputBasePath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii'];
Vvat2.img = VTA_interp; 
ea_write_nii(Vvat2);
ea_autocrop([outputBasePath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii'], '',0,10);
ea_autocrop([outputBasePath, 'efield_model-ossdbs_hemi-', sideLabel, '.nii'], '',0,10);
