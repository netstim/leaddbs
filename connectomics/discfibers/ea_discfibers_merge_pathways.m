function [cfile, map_list, pathway_list] = ea_discfibers_merge_pathways(obj)

% merges pathways from different .mat files and stores them in the
% LeadGroup folder as "merged_pathways.mat"
% also returns global indices of the first fibers in pathways and the
% corresponding list of pathways' names

%myDir = [ea_getconnectomebase('dMRI'), obj.connectome]; % temp
myDir = [ea_getconnectomebase('dMRI_multitract'), obj.connectome];
myFiles = dir(fullfile(myDir,'*.mat')); %gets all mat files in struct

glob_index = 1;
map_list = []; % contains global indices of the first fibers in pathways
pathway_list = {}; % ordered list of pathways' names

C = cell(1,numel(myFiles));
C_idx = cell(1,numel(myFiles));

disp('Merging different pathways ...')

for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  fullFileName = fullfile(myFiles(k).folder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  
  map_list = [map_list, glob_index];
  pathway_list{k} = baseFileName;
  
  fiber_file = load(fullFileName);
  num_of_fibers = length(fiber_file.idx);
  fiber_file.fibers(:,4) = fiber_file.fibers(:,4) + glob_index - 1;
  
  C{k} = fiber_file.fibers;
  C_idx{k} = fiber_file.idx;
  
  glob_index = glob_index + num_of_fibers;
  
end

ftr = fiber_file; % just initialization
% merge cell contents along axis 0
ftr.fibers = cat(1, C{:});
ftr.idx = cat(1, C_idx{:});

if obj.M.ui.detached
    pthprefix = [fileparts(obj.leadgroup),filesep];
else
    pthprefix = '';
end

% store the merged pathways in the leadgroup folder for now
[filepath,name,ext] = fileparts(obj.leadgroup);
cfile = [filepath,filesep,'merged_pathways.mat'];
save(cfile, '-struct', 'ftr');

if obj.statmetric ~= 3
    return
end

% now iterate over fiberActivation.._...mat and merge them
% also adds 0 activation for those filtered out by Kuncel-VTA

numPatient = length(obj.allpatients);
pamlist = cell(numPatient,2);   % no mirroring

disp('Merging fiberActivation files ...')

C_fibState = cell(1,numel(myFiles));
C_fibState_idx = cell(1,numel(myFiles));

for sub=1:numPatient
    for side = 1:2 % hardcoded
        for k=1:length(myFiles)
            
            if side == 1
                side_name = 'right';
            else
                side_name = 'left';
            end
            
            fiberActivation_file = ['fiberActivation_',side_name, '_', myFiles(k).name];
            
            pam_file = [pthprefix, obj.allpatients{sub},filesep, 'stimulations',filesep,...
                ea_nt(0), 'gs_',obj.M.guid,filesep, fiberActivation_file];
            
            
            % we need to add filtered out fibers as not activated
            fib_state_raw = load(char(pam_file));
            
            
            %load(cfile, 'fibers', 'idx');
            total_fibers = length(C_idx{k});
            fib_state = zeros(total_fibers,1);
            
            C_fibState{k} = C{k};
            
            last_loc_i = 1;
            sub_i = 1;
            last_glob = 1;
            for fib_i = 1:total_fibers
                % if the fiber was processed in OSS-DBS, check the status
                if fib_state_raw.fibers(last_loc_i,4) == fib_i
                    fib_state(fib_i) = fib_state_raw.fibers(last_loc_i,5);
                    last_loc_i = fib_state_raw.idx(sub_i)+last_loc_i;
                    sub_i = sub_i + 1;
                    
                else
                    fib_state(fib_i) = 0;  % the fiber was pre-filtered out with Kuncel-VTA
                end
                
                end_index = last_glob+C_idx{k}(fib_i)-1;
                
                C_fibState{k}(last_glob:end_index,5) = fib_state(fib_i);
                
                last_glob = end_index + 1;
                
            end
            
            %C_fibState{k} = [C{k},fib_state];
            C_fibState_idx{k} = C_idx{k};
            
        end
        
        ftr2 = fib_state_raw; % just initialization
        % merge cell contents along axis 0
        ftr2.fibers = cat(1, C_fibState{:});
        ftr2.idx = cat(1, C_fibState_idx{:});
        
        % store as fiberActivation_side.mat in the corresp. stim folder
        [filepath,name,ext] = fileparts(pam_file);
        fiberActivation_merged = [filepath,filesep,'fiberActivation_',side_name,'.mat'];
        save(fiberActivation_merged, '-struct', 'ftr2');
        
    end
end
end