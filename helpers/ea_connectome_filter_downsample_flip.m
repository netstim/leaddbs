function ea_connectome_filter_downsample_flip(connectome_path, ROI_file, factor, fixed_N, N_threshold, flip)
% This function allows to prepare large connnectomes for simulations:
% 1) Optionally, filter OUT fibers that do NOT pass the ROI
% 2) Optionally, downsample the connectome either by factor (integer,
% default option) or to a fixed number. 
% 3) Optionally, flip fibers to another hemisphere (then stored in dMRI_MultiTract)

% By K.Butenko

arguments
    connectome_path                 % path to the folder containing the connectome files
    ROI_file = false                % path to a nifti file containing a binary ROI that has to be intersected by the fibers to be preserved. Use a right hemisphere ROI, when flipping is used
    factor    {mustBeNumeric} = 5   % downsamling factor
    fixed_N    {mustBeNumeric} = 0  % if not 0, downsample to this value
    N_threshold  {mustBeNumeric} = 0  % do not downsample pathways with fibers less than the threshold
    flip   {mustBeNumericOrLogical} = false  % flip fibers to another hemisphere 
end

% check if the downsampling parameter was provided
% and create an output folder

% if flip is active, then save to dMRI_MultiTract
C = strsplit(connectome_path,filesep);
connectome_type = C{end-1};
if flip && strcmp(connectome_type,'dMRI')
    new_connectome_path = strrep(connectome_path,'dMRI','dMRI_MultiTract');
else
    new_connectome_path = connectome_path;
end

if factor ~= 0
    new_connectome_path = [new_connectome_path,'FilteredByROIDownsampledBy',char(string(factor))];
elseif fixed_N ~= 0
    new_connectome_path = [new_connectome_path,'FilteredByROIDownsampledTo',char(string(fixed_N))];
else
    new_connectome_path = [new_connectome_path,'FilteredByROI'];
end
mkdir(new_connectome_path)

if ROI_file
    % load ROI to filter fibers
    ROI = ea_load_nii(ROI_file);
end

%gets all mat files in struct
myFiles = dir(fullfile(connectome_path,'*.mat'));
% remove adjacency matrix and dataset_info if present
myFiles = myFiles(~endsWith({myFiles.name}, '_ADJ.mat'));
myFiles = myFiles(~endsWith({myFiles.name}, '_info.mat'));

for k = 1:length(myFiles)
    baseFileName = myFiles(k).name;
    pathway_file = fullfile(myFiles(k).folder, baseFileName);

    new_pathway = [new_connectome_path,filesep,myFiles(k).name];
    %fprintf(1, 'Now reading %s\n', fullFileName);

    ftr_full = load(pathway_file);

    if ROI_file
        % Trim connectome fibers by ROI
        ROI_Ind = find(abs(ROI.img(:))>0.5);   % ROI is assumed to be binary
    
        % Trim connectome fibers
        [xvox, yvox, zvox] = ind2sub(size(ROI.img), ROI_Ind);
        ROImm = ea_vox2mm([xvox, yvox, zvox], ROI.mat);
        filter = all(ftr_full.fibers(:,1:3)>=min(ROImm),2) & all(ftr_full.fibers(:,1:3)<=max(ROImm), 2);
    
        % discard the pathway if completely unconnected
        if ~any(filter)
            continue;
        end
    
        trimmedFiber = ftr_full.fibers(filter,:);
    
        [trimmedFiberInd, ~, trimmedFiberID] = unique(trimmedFiber(:,4), 'stable');
        fibVoxInd = splitapply(@(fib) {ea_mm2uniqueVoxInd(fib, ROI)}, trimmedFiber(:,1:3), trimmedFiberID);
        
        % Remove outliers
        fibVoxInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
        trimmedFiberInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
      
        trimmedIdx = ftr_full.idx(trimmedFiberInd,:);
        % restore complete trimmed fibers
        trimmedFiber = ftr_full.fibers(ismember(ftr_full.fibers(:,4), trimmedFiberInd), :);
    else
        trimmedFiber = ftr_full.fibers;
        trimmedIdx = ftr_full.idx;
    end

%     segm_counter = 1;
%     for fib_i = 1:length(trimmedFiberInd)
%         fiber_to_copy = ftr_full.fibers(:,4) == trimmedFiberInd(fib_i);
%         trimmedFiber(segm_counter:segm_counter + trimmedIdx(fib_i)-1,:) = ftr_full.fibers(fiber_to_copy,:);
%         segm_counter = segm_counter + trimmedIdx(fib_i);
%     end


    % check if downsampling needed
    if (factor ~= 0 || fixed_N ~= 0) && length(trimmedIdx) > N_threshold
        
        if factor ~= 0
            New_N = round(length(trimmedIdx) / factor);
        else
            New_N = fixed_N;
        end
        
        indices_picked = randperm(length(trimmedIdx),New_N)';
        indices_picked = sort(indices_picked);

        % check how well this works for large dMRI connectomes
        ftr.idx = zeros(length(indices_picked),1);
        ftr.fibers = zeros(sum(trimmedIdx(indices_picked)),4);

        segm_counter = 1;
        for fib_i = 1:length(indices_picked)

            if indices_picked(fib_i) == 1
                start_point = 1;
            else
                start_point = sum(trimmedIdx(1:indices_picked(fib_i)-1)) + 1;
            end
            fib = trimmedFiber(start_point:start_point+trimmedIdx(indices_picked(fib_i))-1,1:3);

            ftr.idx(fib_i) = size(fib,1);
            ftr.fibers(segm_counter:segm_counter + ftr.idx(fib_i) - 1,1:3) = fib;
            ftr.fibers(segm_counter:segm_counter + ftr.idx(fib_i) - 1,4) = fib_i;
            segm_counter = segm_counter + ftr.idx(fib_i);
        end

    else   % no downsampling, just copy the file, but adjust indexing after ROI filtering
        ftr.fibers = zeros(sum(trimmedIdx),4);
        ftr.idx = trimmedIdx;

        orig_indices = unique(trimmedFiber(:,4));
        for inx = 1:length(orig_indices)
            inx_to_change = trimmedFiber(:,4) == orig_indices(inx);
            trimmedFiber(inx_to_change,4) = inx;
        end

        ftr.fibers = trimmedFiber;

    end

    % save the result
    ftr.fourindex = 1;
    ftr.ea_fibformat = '1.0';
    save(new_pathway, '-struct', 'ftr');

    if flip

        ftr.mirrored = 1;  % this flag designates that the connectome is fiberwise mirrorred
        save(new_pathway, '-struct', 'ftr');

        % now flip and store (we assume that the flip is from the right)
        fibers_coords = ftr.fibers(:,1:3);
        fibers_coords_lr = ea_flip_lr_nonlinear(fibers_coords);
        ftr.fibers(:,1:3) = fibers_coords_lr;
        
        % add left hemisphere suffix
        new_pathway = [new_connectome_path,filesep,myFiles(k).name(1:end-4), '_flipped', '.mat'];
        save(new_pathway, '-struct', 'ftr');
        
    end

end

% % add mirror flag to already processed connectomes
% %gets all mat files in struct
% myFiles = dir(fullfile(connectome_path,'*.mat'));
% % remove adjacency matrix if present
% myFiles = myFiles(~endsWith({myFiles.name}, '_ADJ.mat'));
% 
% for k = 1:length(myFiles)
%     baseFileName = myFiles(k).name;
%     pathway_file = fullfile(myFiles(k).folder, baseFileName);
%     ftr = load(pathway_file);
% 
%     ftr.mirrored = 1;  % this flag designates that the connectome is fiberwise mirrorred
%     save(pathway_file, '-struct', 'ftr');
% end
