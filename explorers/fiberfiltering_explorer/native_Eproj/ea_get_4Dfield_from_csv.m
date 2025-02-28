function ea_get_4Dfield_from_csv(field_in_csv, file2save)
% Get nifti file with X,Y,Z vector components of the electric field.
% By Butenko, konstantinmgtu@gmail.com

arguments
    field_in_csv     % where the csv file with field components is stored
    file2save        % 4D nifti filename
end

%ea_dispt('Creating nifti header for export...');
% create nifti
[~, ~, endian] = computer;
switch endian
    case 'L'
        endian = 0;
    case 'B'
        endian = 1;
end

% split to corrdinates and field
Field_array = table2array(readtable(field_in_csv));
Field_coords = Field_array(:,2:4);
n_points = 71;
Field_vals = zeros(n_points,n_points,n_points,3);

% resample the field to a regular grid
for d=1:3  % iterate over field components
    %ea_dispt('Converting to equispaced image data...');
    F = scatteredInterpolant(Field_coords(:,1),Field_coords(:,2),Field_coords(:,3),Field_array(:,4+d),'linear','none');
    gv=cell(3,1); spacing=zeros(3,1);
    
    for axis = 1:3
        gv{axis}=linspace(min(round(Field_coords(:,axis))),max(round(Field_coords(:,axis))),n_points);
        spacing(axis)=abs(gv{axis}(1)-gv{axis}(2)); 
    end
    Field_vals(:,:,:,d) = F(gv);
end

% cap values
Field_vals(isnan(Field_vals)) = 0;
Field_vals(Field_vals>1000.0) = 1000.0;

% I have no idea what is happening here
chun1=randperm(n_points); chun2=randperm(n_points); chun3=randperm(n_points); 
Field_interp.mat=mldivide([(chun1);(chun2);(chun3);ones(1,n_points(1))]',[gv{1}(chun1);gv{2}(chun2);gv{3}(chun3);ones(1,n_points)]')';

% store as a 4-D nifti
Field_interp.img = Field_vals;
Field_interp.dim=[n_points,n_points,n_points];
Field_interp.dt = [4, endian];
Field_interp.n=[1 1];
Field_interp.descrip='oss-dbs-v2 - Field_ref';
Field_interp.fname = file2save;
till_save_nii([Field_interp,Field_interp,Field_interp],Field_interp.img);
