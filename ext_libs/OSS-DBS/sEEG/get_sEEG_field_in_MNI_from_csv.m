function get_sEEG_field_in_MNI_from_csv(field_in_csv, file2save, phi_field, anchor_img, transform, zipIt)
% Get nifti of field
% By Butenko, konstantinmgtu@gmail.com

arguments
    field_in_csv     % where the csv file with field components is stored
    file2save        % path for niftis to save
    phi_field        % if true, potential is exported
    anchor_img       % path to the anchor image
    transform        % from native to MNI, you need the opposite warp, i.e. from-MNI152NLin2009bAsym_to-anchorNative
    zipIt {mustBeNumericOrLogical} = true % gzip if true
end

% split to coordinates and field
Field_array = table2array(readtable(field_in_csv));
Field_coords = Field_array(:,2:4);

% convert to native voxel space (will be as floating numbers)
Field_vox_native = ea_mm2vox(Field_coords, anchor_img)';

ea_dispt('Transforming spatial coordinates to MNI space...');
Field_coords_MNI  = ea_map_coords(Field_vox_native, ...
    anchor_img, ...
    transform, ...
    [ea_space,'t1.nii'], 'ANTS')';

%ea_dispt('Creating nifti header for export...');
% create nifti
[~, ~, endian] = computer;
switch endian
    case 'L'
        endian = 0;
    case 'B'
        endian = 1;
end

if phi_field
    F = scatteredInterpolant(Field_coords_MNI(:,1),Field_coords_MNI(:,2),Field_coords_MNI(:,3),Field_array(:,5),'linear','none');
else
    F = scatteredInterpolant(Field_coords_MNI(:,1),Field_coords_MNI(:,2),Field_coords_MNI(:,3),Field_array(:,8)*1000.0,'linear','none');
end
gv=cell(3,1);

% hardwired N of points, if changed, also change Lattice shape in lead_settings.py
n_points = 71;
for axis = 1:3
    gv{axis}=linspace(min(round(Field_coords_MNI(:,axis))),max(round(Field_coords_MNI(:,axis))),n_points);
end
chun1=randperm(n_points); chun2=randperm(n_points); chun3=randperm(n_points); 
ROI.mat=mldivide([(chun1);(chun2);(chun3);ones(1,n_points(1))]',[gv{1}(chun1);gv{2}(chun2);gv{3}(chun3);ones(1,n_points)]')';
ROI.dim=[n_points,n_points,n_points];

eeg = F(gv);
eeg(isnan(eeg))=0;

%% Export stuff
ea_dispt('Writing files...');

ROI.fname=file2save;
if phi_field
    ROI.descrip='OSS Recording Field';
else
    ROI.descrip='OSS E-Field';
end
ROI.img=eeg; %permute(eeg,[2,1,3]);
ROI.pinfo = [1;0;352];
ROI.dt = [64, endian];
ROI.n=[1 1];
ea_write_nii(ROI);
if zipIt
    gzip(ROI.fname)
    ea_delete(ROI.fname)
end
