function ea_downsample_connectome(connectome_path, factor, fixed_N, N_threshold)

% uniformly downsample connectomes either by factor (integer, default option) or to a fixed
% number (set factor = 0 and fixed_N to the desired number)

% connectome name - could be dMRI or dMRI_MultiTract connectome
% factor - downsamling factor
% N_threshold - only downsample pathways with N > N_threshold

% downsample all pathways if no threshold passed
if ~exist('N_threshold','var')
    N_threshold = 0;
end

% check that the downsampling parameter was provided
% and create an output folder
if factor == 0 && fixed_N == 0
    ea_warndlg("Either factor or fixed_N should be non-zero")
    return
elseif factor ~= 0
    new_connectome_path = [connectome_path,'_downsampled_by_',char(string(factor))];
else
    new_connectome_path = [connectome_path,'_downsampled_to_',char(string(factor))];
end
mkdir(new_connectome_path)

%gets all mat files in struct
myFiles = dir(fullfile(connectome_path,'*.mat'));

% remove adjacency matrix if present
myFiles = myFiles(~endsWith({myFiles.name}, '_ADJ.mat'));


for k = 1:length(myFiles)
    baseFileName = myFiles(k).name;
    pathway_file = fullfile(myFiles(k).folder, baseFileName);
    new_pathway = [new_connectome_path,filesep,myFiles(k).name];
    %fprintf(1, 'Now reading %s\n', fullFileName);

    ftr_full = load(pathway_file);

    % check if downsampling needed
    if length(ftr_full.idx) > N_threshold
        
        if factor ~= 0
            New_N = round(length(ftr_full.idx) / factor);
        else
            New_N = N_threshold;
        end
        
        indices_picked = randperm(length(ftr_full.idx),New_N)';
        indices_picked = sort(indices_picked);

        % check how well this works for large dMRI connectomes
        ftr.idx = zeros(length(indices_picked),1);
        ftr.fibers = zeros(sum(ftr_full.idx(indices_picked)),4);

        segm_counter = 1;
        for fib_i = 1:length(indices_picked)

            fib = ftr_full.fibers(ftr_full.fibers(:,4) == indices_picked(fib_i),:);
            ftr.idx(fib_i) = size(fib,1);
            ftr.fibers(segm_counter:segm_counter + ftr.idx(fib_i) - 1,:) = fib;
            ftr.fibers(segm_counter:segm_counter + ftr.idx(fib_i) - 1,4) = fib_i;
            segm_counter = segm_counter + ftr.idx(fib_i);

%             idx_fib_i = 0;
%             % instead of searching, just go block-by-block
%             while ftr.fibers(segm_counter, 4) == indices_picked(fib_i)
%                 ftr.fibers(segm_counter, 4) = fib_i;  % re-index
%                 segm_counter = segm_counter + 1;
%                 idx_fib_i = idx_fib_i + 1;
%             end
        end

        ftr.fourindex = 1;
        ftr.ea_fibformat = '1.0';
        save(new_pathway, '-struct', 'ftr');

    else   % no downsampling, just copy the file
        copyfile(pathway_file, new_pathway)
    end
end
