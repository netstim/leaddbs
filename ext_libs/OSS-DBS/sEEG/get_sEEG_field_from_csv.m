function get_sEEG_field_from_csv(field_in_csv, file2save, phi_field, reslice2segmask, segmaskFile)
% Get nifti of OSS electric fields
% By Butenko, konstantinmgtu@gmail.com

arguments
    field_in_csv     % where the csv file with field components is stored
    file2save        % path for niftis to save
    phi_field        % if true, potential is exported
    reslice2segmask  % if true, all VTRs are stored in the same voxel space defined by segmask
    segmaskFile      % path to the nii
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
if phi_field
    F = scatteredInterpolant(Field_coords(:,1),Field_coords(:,2),Field_coords(:,3),Field_array(:,5),'linear','none');
else
    F = scatteredInterpolant(Field_coords(:,1),Field_coords(:,2),Field_coords(:,3),Field_array(:,8)*1000.0,'linear','none');
end
gv=cell(3,1);

if reslice2segmask
    res = 0.5;
    segmaskFileResliced = [segmaskFile(1:end-4),'_resliced.nii'];
    ea_reslice_nii(segmaskFile, segmaskFileResliced, res)
    segmask = ea_load_nii(segmaskFileResliced);
    n_points = round(segmask.dim);  % change this for higher resolution
    first_point = segmask.mat * [0,0,0,1]';
    last_point = segmask.mat * [segmask.dim(1),segmask.dim(2),segmask.dim(3),1]';
    segmask_res = segmask;
    segmask_res.dim = n_points;
    segmask_res.img = zeros(n_points);
    
    for axis = 1:3
        gv{axis}=linspace(first_point(axis,1),last_point(axis,1),n_points(axis));
    end
    ROI = segmask_res;
else
    % hardwired N of points, if changed, also change Lattice shape in lead_settings.py
    n_points = 71;
    for axis = 1:3
        gv{axis}=linspace(min(round(Field_coords(:,axis))),max(round(Field_coords(:,axis))),n_points);
    end
    chun1=randperm(n_points); chun2=randperm(n_points); chun3=randperm(n_points); 
    ROI.mat=mldivide([(chun1);(chun2);(chun3);ones(1,n_points(1))]',[gv{1}(chun1);gv{2}(chun2);gv{3}(chun3);ones(1,n_points)]')';
    ROI.dim=[n_points,n_points,n_points];
end

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
gzip(ROI.fname)
ea_delete(ROI.fname)

