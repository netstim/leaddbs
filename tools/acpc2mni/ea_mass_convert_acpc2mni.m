function ea_mass_convert_acpc2mni(mass_convert_config, input_table_filepath, uipatdir, output_save_filepath)

% read the input ACPC coordinates
A = readtable(input_table_filepath);

acpc_coords = A.InputAC_PCCoordinate;
acpc_coords = split(acpc_coords);
acpc_coords = str2double(acpc_coords);


% negate the acpc coordinate if necessary
if strcmp(mass_convert_config.lr, 'left') % TODO add checks to ensure that this is either set to "left" or "right"
    acpc_coords(:,1) = -1 * acpc_coords(:,1);
end
if strcmp(mass_convert_config.ap, 'posterior')
    acpc_coords(:,2) = -1 * acpc_coords(:,2);
end
if strcmp(mass_convert_config.ab, 'below')
    acpc_coords(:,3) = -1 * acpc_coords(:,3);
end

all_cfgs = acpc_coords; % copy over acpc coordinates
all_cfgs = [all_cfgs zeros(height(A),1)]; % mapmethod
all_cfgs = [all_cfgs 2*ones(height(A),1)]; % acmcpc


mni_coords = zeros(height(A),3);
parfor row = 1:height(A)
    disp(A(row,:))
    cfg = struct()
    cfg.xmm = all_cfgs(row, 1);
    cfg.ymm = all_cfgs(row, 2);
    cfg.zmm = all_cfgs(row, 3);
    cfg.mapmethod = all_cfgs(row, 4);
    cfg.acmcpc = all_cfgs(row, 5);
    fid1=ea_acpc2mni(cfg,uipatdir);
%     averaged_coords = average_mni_coordinate_for_one_patient(fid1); % TODO use mean() function instead of custom function. Reference ea_acpcquery.m
    
    
    % average MNI coordinates for one patient
    sum_array = [0,0,0];
    fid_size = size(fid1);

    for i=1:fid_size(2)
        sum_array = sum_array + fid1(i).WarpedPointMNI;
    end

    averaged_coords = sum_array / fid_size(2);
    
    mni_coords(row,:) = averaged_coords
end
disp(mni_coords)

output.acpc_coords = acpc_coords
output.mni_coords = mni_coords
save(output_save_filepath, 'output')
end