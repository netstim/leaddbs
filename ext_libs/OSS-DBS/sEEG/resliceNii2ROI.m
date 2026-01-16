function resliceNii2ROI(images2reslice,threshold,vox_size)
    % Reslice nifti to the common voxel space
    % By Butenko, konstantinmgtu@gmail.com
    
    arguments
        images2reslice    % cell array of full paths, Nx1 
        threshold        % float, cutoff value to define the bounding box
        vox_size         % vector of floats 1x3, voxel dimensions,
    end

    [~, ~, endian] = computer;
    switch endian
        case 'L'
            endian = 0;
        case 'B'
            endian = 1;
    end

    % --- Step 1: Find the Global Bounding Box (World Space) ---
    fprintf('Calculating common bounding box...\n');
    all_min = [Inf Inf Inf];
    all_max = [-Inf -Inf -Inf];
    
    for i = 1:length(images2reslice)
        V = spm_vol(images2reslice{i});
        img = spm_read_vols(V);
        
        % Find coordinates of voxels exceeding threshold
        idx = find(img < threshold);
        if isempty(idx), continue; end
        
        [r, c, s] = ind2sub(size(img), idx);
        vox_coords = [r c s ones(length(r),1)]';
        
        % Convert voxel indices to world coordinates (mm) using the affine matrix
        world_coords = V.mat * vox_coords;
        world_coords = world_coords(1:3, :)';
        
        % Update the global min/max
        all_min = min([all_min; world_coords]);
        all_max = max([all_max; world_coords]);
    end
    
    % Define the final Bounding Box (BB) in SPM format: [minX minY minZ; maxX maxY maxZ]
    BB = [all_min; all_max];
    
    % --- Step 2: Define the Target Affine (New Header) ---
    % We create a dummy header to define the grid for reslicing
    % Origin is set to the minimum corner of the bounding box
    mat = spm_matrix(all_min(:)') * diag([vox_size 1]) * spm_matrix(-[1 1 1]);
    dim = ceil(diff(BB) ./ vox_size) + 1;
    
    % Create a temporary reference volume structure
    target_V = struct('img',zeros(dim),'mat', mat, 'dim', dim, 'fname', 'temp_ref.nii', 'pinfo', [1 0 0]', 'dt', [4, endian], 'n', [1,1]);
    ea_write_nii(target_V);
    
    % --- Step 3: Reslice All Images to this Space ---
    flags = struct('interp', 4, 'mask', 0, 'mean', 0, 'which', 1, 'prefix', 'resliced_');
    
    for i = 1:length(images2reslice)
        % fprintf('Reslicing %s...\n', subj_RFs{i});
        % % We pass [target_V; source_V] to spm_reslice
        % source_V = spm_vol(subj_RFs{i});
        % spm_reslice({target_V, source_V}, flags);
        ea_conformspaceto(target_V.fname,images2reslice{i})

        nii_rescliced = ea_load_nii(images2reslice{i});
        nii_rescliced.img(isnan(nii_rescliced.img)) = 0.0;
        nii_rescliced.fname = [nii_rescliced.fname,'.gz'];
        ea_write_nii(nii_rescliced);
        ea_delete(images2reslice{i})

    end
    %fprintf('All images resliced to the same voxel space.\n'); 
end