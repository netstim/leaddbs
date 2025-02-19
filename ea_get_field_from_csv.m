function ea_get_field_from_csv(ref_image, Field_array_file, Activation_threshold_VTA, sideLabel, outputBasePath, source_index)

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
Field_array = table2array(readtable(Field_array_file, ReadVariableNames=false));
Field_coords = Field_array(:,2:4);
Field_vals = Field_array(:,8) * 1000.0;  % others are the components

% just interpolate the magnitude, the vector field can be confusing
ea_dispt('Converting to equispaced image data...');
F = scatteredInterpolant(Field_coords(:,1),Field_coords(:,2),Field_coords(:,3),Field_vals,'linear','none');
gv=cell(3,1); spacing=zeros(3,1);

% hardwired N of points, if changed, also change Lattice shape in lead_settings.py
n_points = 71;
for axis = 1:3
    %n_points(axis) = (max(round(Field_coords_MNI(:,axis))) - min(round(Field_coords_MNI(:,axis)))) / template.voxsize(axis);
    gv{axis}=linspace(min(round(Field_coords(:,axis))),max(round(Field_coords(:,axis))),n_points);
    spacing(axis)=abs(gv{axis}(1)-gv{axis}(2)); 
end

% I have no idea what is happening here
chun1=randperm(n_points); chun2=randperm(n_points); chun3=randperm(n_points); 
Field_interp.mat=mldivide([(chun1);(chun2);(chun3);ones(1,n_points(1))]',[gv{1}(chun1);gv{2}(chun2);gv{3}(chun3);ones(1,n_points)]')';

Field_interp.dim=[n_points,n_points,n_points];
Field_interp.dt = [4, endian];
Field_interp.n=[1 1];
Field_interp.descrip='oss-dbs-v2 - Field_ref';

[filepath,filename,~] = fileparts(Field_array_file);

Field_interp.img = F(gv);
Field_interp.img(isnan(Field_interp.img)) = 0;
Field_interp.img(Field_interp.img>1000.0) = 1000.0;
if source_index == 5  % no source indexing
    Field_interp.fname = [outputBasePath, 'efield_model-ossdbs_hemi-', sideLabel, '.nii'];
else
    Field_interp.fname = [outputBasePath, 'efield_model-ossdbs_hemi-', sideLabel,'_S',num2str(source_index), '.nii'];
end
ea_write_nii(Field_interp);
%ea_autocrop([outputBasePath, 'efield_model-ossdbs_hemi-', sideLabel, '.nii'], margin=10);

% also create VATs directly
VTA_interp = Field_interp.img >= (Activation_threshold_VTA);
Vvat2 = Field_interp;
Vvat2.descrip='oss-dbs-v2 - VAT_ref';
if source_index == 5  % no source indexing
    Vvat2.fname = [outputBasePath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii'];
else
    Vvat2.fname = [outputBasePath, 'binary_model-ossdbs_hemi-', sideLabel,'_S',num2str(source_index), '.nii'];
end
Vvat2.pinfo = [1;0;352];
Vvat2.dt = [2, endian];
Vvat2.img = VTA_interp; 
ea_write_nii(Vvat2);
%ea_autocrop([outputBasePath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii'],, margin=10);
%ea_autocrop([outputBasePath, 'efield_model-ossdbs_hemi-', sideLabel, '.nii'],, margin=10);
