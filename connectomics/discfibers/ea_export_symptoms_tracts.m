function ea_export_symptoms_tracts(obj, symptomTractsfile)

% this function creates symptom-specific tracts from fiber filtering
% results. Those will be separated based on the spatial distribution and
% the distribution of vals (fiber statistic, e.g.  t- or r-values)

% IMPORTANT: use either positive or negative tracts
% Negative tracts can be used for soft side-effects (their weights will be abs(norm((vals)))

% unlike export fibscore model, we do not need to map to global space,
% but only have the global indices

% we want to define the model outside of cross-validation
if ~exist('patsel','var') % patsel can be supplied directly (in this case, obj.patientselection is ignored), e.g. for cross-validations.
    patientsel = obj.patientselection;
end

% get fiber model (vals) and corresponding indices (usedidx)
if obj.cvlivevisualize
    [vals,fibcell,usedidx] = ea_discfibers_calcstats(obj, patientsel);
    obj.draw(vals,fibcell,usedidx)
    drawnow;
else
    [vals,~,usedidx] = ea_discfibers_calcstats(obj, patientsel);
end

%% parameters

% threshold vicinity mask
vcnty_thr = 15; % N.B: units are in mm (15 mm works well)

% Gaussian Kernel Size applied to fibers when computing spatial correlation
sig = 1.5; 

% threshold for spatial hierarchical clustering (the larger the value, the more groups you get)
prox_cluster_threshold = 0.17; % not yet optimal

var_threshold = 0.01; % acceptable variance of val within one spatial group
val_metric_coef = 0.9; % accepted fraction of otsu metric for val variance (1.0 is defined by the maximum number of bins)

% should be passed to the function
symptoms_list = {'Brady', 'Rigidity', 'Axial', 'Tremor'};

% we extract each voter and side separately
for voter=1:size(vals,1)  % I would restrict to one voter for now
    for side = 1:size(vals,2)

        % target coordinates that define the center of masked brain region where
        % spatial correlation will be computed
        if side == 1
            prefix_side = '_rh';
            trgt_coor = [11.7631,-13.9674,-8.85935]; % target coordinate (e.g STN), right and left
        else
            prefix_side = '_lh';
            trgt_coor = [-11.7631,-13.9674,-8.85935];
        end

        % load cfile
        [filepath,~,~] = fileparts(obj.leadgroup);
        if obj.multi_pathways == 1
            cfile = [filepath,filesep,'merged_pathways.mat'];
        else
            cfile = [ea_getconnectomebase('dMRI'), obj.connectome, filesep, 'data.mat'];
        end
        load(cfile, 'fibers', 'idx');

        % we can check and remove fibers that are too far away from the
        % target (but track the indices!)
        
        if obj.connectivity_type == 2
            gl_indices = obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{side}(usedidx{voter,side});
        else
            gl_indices = obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{side}(usedidx{voter,side});
        end

        % normalized to 0-1 if necessary (apply abs if negative)
        if any(min(vals{voter,side}) < 0) && any(max(vals{voter,side}) > 0)
            disp("Choose only positive or negative tracts!")
            return;
        elseif min(vals{voter,side}) < 0
            vals{voter,side} = abs(vals{voter,side});
            % those are for soft-threshold side-effects
            negative_vals = 1;
        else
            negative_vals = 0;
        end

        if max(vals{voter,side}) > 1.0
            vals{voter,side} = vals{voter,side} / max(vals{voter,side});
        end


        % think about some standard way of symptom naming
        if size(vals,1) > 1
            % think about a function here for multitracts,
            symptomName = [obj.subscore.labels{voter,1},prefix_side];
        else
            symptomName = [obj.responsevarlabel,prefix_side];
        end
    
        % change - to _
        symptomName = strrep(symptomName,'-','_');

        % prepare output folders for pathways
        symptomNameFolder = [filepath,filesep,symptomName(1:end-3)];

        % store soft side-effect pathways separately
        if negative_vals == 1
            symptomNameFolder = [symptomNameFolder,'_SE'];
        end

        if side == 1 
            if isfolder(symptomNameFolder)
                rmdir(symptomNameFolder, 's')
            end
            mkdir(symptomNameFolder)
        end        


        % create fibers and idx only from remaining fibers (with vals)
        % but remember to use gl_indices
        fibers_with_vals = [];
        fibers_orig_ind = [];
        idx_with_vals = idx(gl_indices);
        fibers_all = fibers;
        idx_all = idx;

        for i = 1:length(gl_indices)
            fibers_with_vals = cat(1,fibers_with_vals,fibers_all(fibers_all(:,4) == gl_indices(i),:));
            fibers_with_vals(fibers_with_vals(:,4) == gl_indices(i),4) = i;
            fibers_orig_ind = cat(1,fibers_orig_ind,fibers_all(fibers_all(:,4) == gl_indices(i),4));
        end

        % Get spatial correlation matrix for fibers after Gaussian filter
        % should be accelerated with parfor
        spatio_corr_mat = ea_get_fiber_spatial_corr(fibers_with_vals,trgt_coor,sig,vcnty_thr);

        % leave only upper triangle values (the matrix is symmetric)
        PRX_M = triu(spatio_corr_mat,1)';

        %% hierarchical clustering of the correlation matrix
        
        % you can check mean(PRX_M) to decide on the clustering threhsold

        % and convert to a vector (as pdist)
        dissimilarity = 1 - PRX_M(find(PRX_M))';

        % remember that 0.4 corresponds to corr of 0.6!
        cutoff = 1 - prox_cluster_threshold; 

        % perform linkage clustering
        Z = linkage(dissimilarity);
        groups = cluster(Z,'cutoff',cutoff,'criterion','distance');

        % we can visualize groups using OSS-DBS output format
%         for i = 1:length(idx_with_vals)
%             fibers_with_vals(fibers_with_vals(:,4) == i,5) = groups(i)-4;
%         end
%         origNum = length(idx);
%         connectome_name = 'placeholder';

        % rename just to be able to load later
%         fibers = fibers_with_vals; % IMPORTANT: remember that 4th column are local indices here!
%         idx = idx_with_vals;
        % then save to mat origNum, fibres, connectome_name and idx

        % decrease number of value bins if a lot of spatial groups
        if groups > 10
            disp("Warning: number of spatially distributed pathways above 10")
            disp("Number of val bins per pathway will be reduced to 5")
            max_N_pathways = 5; % per spatially defined pathway
        else
            max_N_pathways = 10;
        end


        % Now bin within blocks based on val variance 
        pathway_index_list = cell(1,length(unique(groups)));
        pathway_mean_vals = cell(1,length(unique(groups)));
        % sum val over the hemisphere to weight predictions across
        % hemispheres
        pathway_sum_vals = cell(1,length(unique(groups)));

        disp("Number of spatial groups "+string(length(unique(groups))))

        % we remove pathways with less than 10 fibers, because they are 
        % not suitable for the analysis 
        drop_pathway = zeros(length(unique(groups)), 1);

        for group_i = 1:length(unique(groups))

            if length(vals{voter,side}(groups == group_i)) < 10
                drop_pathway(group_i) = 1;
                continue
            end

            % do not attempt to val bin one fiber
            if length(vals{voter,side}(groups == group_i)) > 1
    
                if max_N_pathways > length(find(groups == group_i))
                    max_N_pathways_bin = length(find(groups == group_i));
                else
                    max_N_pathways_bin = max_N_pathways;
                end
    
                %metric_vals = zeros(max_N_pathways_bin-1,1);

                % start with the max number of bins, and check how worse
                % other perform based on metric
                for class_number = max_N_pathways_bin:-1:1
                    
                    input = vals{voter,side}(groups == group_i);  
                    idx_group = gl_indices(groups == group_i);
                        
                    % check variance within the group
                    % do not split if too small
                    % also fixes the issue of one fiber
                    if var(input) < var_threshold
                        thresh = mean(input);
                        metric = -1;
                        break
                    end

                    % Otsu's variance based method
                    [thresh, metric] = multithresh(input,class_number);
                    if metric == 0 && length(input) > 1
                        % did not work, do not split
                        thresh = mean(input);
                        metric = -1;
                        break
                    end

                    %metric_vals(class_number-1) = metric;
                    if class_number == max_N_pathways_bin
                        metric10 = metric; % highest metric
                    else
                        if metric < val_metric_coef * metric10 % check this condition empirically
                            fprintf("N classes: %d \n", class_number + 1 )
                            disp(metric)
                            disp(metric10)                            
                            % recompute thresholds
                            [thresh, metric] = multithresh(input,class_number + 1);
                            break
                        end
                    end
                end
    
    
                if metric == -1 % only one val per spatial pathway
                    pathway_index_list{1,group_i}{1,1} = idx_group;
                    pathway_mean_vals{1,group_i}(1) = mean(input); 
                    pathway_sum_vals{1,group_i}(1) = sum(input);
                else
                    pathway_index_list{1,group_i} = cell(1,length(thresh) + 1);
                    pathway_mean_vals{1,group_i} = zeros(1,length(thresh) + 1);
                    pathway_sum_vals{1,group_i} = zeros(1,length(thresh) + 1);

                    % group to bins based on thresholds
                    for thresh_i = 1:length(thresh) + 1
        
                        if thresh_i == 1
        
                            pathway_mean_vals{1,group_i}(thresh_i) = mean(input(input < thresh(thresh_i)));
                            pathway_sum_vals{1,group_i}(thresh_i) = sum(input(input < thresh(thresh_i)));
                            pathway_index_list{1,group_i}{thresh_i} = idx_group(input < thresh(thresh_i));
                                
%                             % alternatively, drop pathways with vals < 0.5 and N fibers < 10
%                             if pathway_mean_vals{1,group_i}(thresh_i) < 0.5 && length(find(idx_group(input < thresh(thresh_i)))) < 10
%                                 pathway_index_list{1,group_i}{thresh_i} = 0;
%                             else
%                                 pathway_index_list{1,group_i}{thresh_i} = idx_group(input < thresh(thresh_i));
%                             end
        
                        elseif thresh_i == length(thresh) + 1
        
                            pathway_mean_vals{1,group_i}(thresh_i) = mean(input(input >= thresh(thresh_i-1)));
                            pathway_sum_vals{1,group_i}(thresh_i) = sum(input(input >= thresh(thresh_i-1)));
                            pathway_index_list{1,group_i}{thresh_i} = idx_group(input >= thresh(thresh_i-1));

        %                     % alternatively, drop pathways with vals < 0.5 and N fibers < 10
        %                     if pathway_mean_vals{1,group_i}(thresh_i) < 0.5 && length(pathway_index_list{1,group_i}{thresh_i}) < 10
        %                         pathway_index_list{1,group_i}{thresh_i} = 0;
        %                     else
        %                         pathway_index_list{1,group_i}{thresh_i} = gl_indices(find(input > thresh(thresh_i-1)));
        %                     end
                        else 
                            pathway_mean_vals{1,group_i}(thresh_i) = mean(input(input >= thresh(thresh_i-1) & input < thresh(thresh_i)));
                            pathway_sum_vals{1,group_i}(thresh_i) = sum(input(input >= thresh(thresh_i-1) & input < thresh(thresh_i)));
                            pathway_index_list{1,group_i}{thresh_i} = idx_group(input >= thresh(thresh_i-1) & input < thresh(thresh_i));

%                             % alternatively, drop pathways with vals < 0.5 and N fibers < 10
%                             if pathway_mean_vals{1,group_i}(thresh_i) < 0.5 && length(find(input >= thresh(thresh_i-1) & input < thresh(thresh_i))) < 10
%                                 pathway_index_list{1,group_i}{thresh_i} = 0;
%                             else
%                                 pathway_index_list{1,group_i}{thresh_i} = idx_group(input >= thresh(thresh_i-1) & input < thresh(thresh_i));
%                             end
        
                        end
                    end
                end
            else
                % if one fiber per spatial pathway
                pathway_mean_vals{1,group_i}(1) = vals{voter,side}(groups == group_i);
                pathway_sum_vals{1,group_i}(1) = vals{voter,side}(groups == group_i);
                pathway_index_list{1,group_i}{1,1} = gl_indices(groups == group_i);
            end

            % we might consider more filtering here

            % save as a pathway
            %fibers, idx, fourindex, ea_fibformat and fibers_glob_index
            ea_fibformat = '1.0';
            fourindex = 1;

            % iterate over spatial groups
            for i= 1:length(pathway_index_list{1,group_i})
                fibers_pathway = [];
                fibers_glob_ind = [];

                if length(pathway_index_list{1,group_i}{1,i}) == 1
                    fibers_pathway = fibers_all(fibers_all(:,4) == pathway_index_list{1,group_i}{1,i},:);
                    fibers_pathway(:,4) = 1;
                    fibers_glob_ind = pathway_index_list{1,group_i}{1};
                    idx_pathway = size(fibers_pathway,1);
                else
                    idx_pathway = zeros(length(pathway_index_list{1,group_i}{1,i}),1);
                    % iterative over val bins
                    for j = 1:length(pathway_index_list{1,group_i}{1,i})
                        fibers_pathway = cat(1,fibers_pathway,fibers_all(fibers_all(:,4) == pathway_index_list{1,group_i}{1,i}(j),:));
                        fibers_pathway(fibers_pathway(:,4) == pathway_index_list{1,group_i}{1,i}(j),4) = j;
                        fibers_glob_ind = cat(1,fibers_glob_ind,fibers_all(fibers_all(:,4) == pathway_index_list{1,group_i}{1,i}(j),4));
                        idx_pathway(j) = size(fibers_all(fibers_all(:,4) == pathway_index_list{1,group_i}{1,i}(j),4),1);
                    end
                end
                ftr.fibers = fibers_pathway;
                ftr.idx = idx_pathway;
                ftr.ea_fibformat = ea_fibformat;
                ftr.fourindex = fourindex;
                % store metadata
                ftr.glob_ind = fibers_glob_ind;
                ftr.orig_connectome = obj.connectome;

                if negative_vals == 0
                    pathway_name = char(compose('%s_%s_group_%d_val_%.3f',symptoms_list{voter},prefix_side,group_i,pathway_mean_vals{1,group_i}(i)));
                else
                    pathway_name = char(compose('SE_%s_%s_group_%d_val_%.3f',symptoms_list{voter},prefix_side,group_i,pathway_mean_vals{1,group_i}(i)));
                end
                
                if contains(pathway_name,'1.0')
                    pathway_name = erase(pathway_name,'.0'); % remove digits after comma
                else   
                    pathway_name = strrep(pathway_name,'0.','0');
                end

                % you should also add the symptom name
                save([symptomNameFolder,filesep,pathway_name],'-struct','ftr')
            end

        end

        % in FF, it is natural that the target activation for positive
        % tracts is 1.0 and for negative (soft-side effect) it is 0.0

        % weights are determined by pathway_mean_vals

        % we iterate over sides and voters and add symptoms to the
        % file
        % if the tract is already in the dict, it will be overwritten

        % load previous if available
        if isfile(symptomTractsfile)
            jsonText = fileread(symptomTractsfile);
            % Convert JSON formatted text to MATLAB data types 
            jsonDict = jsondecode(jsonText); 
        end

        for group_i = 1:size(pathway_mean_vals,2)

            % skip dropped
            if drop_pathway(group_i) == 1
                continue
            end

            for bin_i = 1:length(pathway_mean_vals{1,group_i})

                % remove period delimiter
                if negative_vals == 0
                    pathway_name = char(compose('%s_%s_group_%d_val_%.3f',symptoms_list{voter},prefix_side,group_i,pathway_mean_vals{1,group_i}(bin_i)));
                else
                    pathway_name = char(compose('SE_%s_%s_group_%d_val_%.3f',symptoms_list{voter},prefix_side,group_i,pathway_mean_vals{1,group_i}(bin_i)));
                end

                if contains(pathway_name,'1.0')
                    pathway_name = erase(pathway_name,'.0'); % remove digits after comma
                else   
                    pathway_name = strrep(pathway_name,'0.','0');
                end

                if negative_vals == 0
                    jsonDict.profile_dict.(genvarname(symptomName)).(genvarname(pathway_name)) = [1.0, pathway_mean_vals{1,group_i}(bin_i), pathway_sum_vals{1,group_i}(bin_i)];
                else
                    jsonDict.Soft_SE_dict.(genvarname(['SE_',symptomName])).(genvarname(pathway_name)) = [0.0, pathway_mean_vals{1,group_i}(bin_i), pathway_sum_vals{1,group_i}(bin_i)];
                end
            end
        end

        if all(drop_pathway == 1)
            ea_warndlg("No pathways were exported for the symptom")
        else
            % save to the same file
            jsonText2 = jsonencode(jsonDict);
            fid = fopen(symptomTractsfile, 'w');
            fprintf(fid, '%s', jsonText2)
            fclose(fid);
        end

    end
end        

if sum(drop_pathway) > 0
    ea_warndlg("Pathways with less than 10 fibers were removed, check the terminal")
    disp("Number of small pathways")
    disp(sum(drop_pathway))
end

% check number or repetitions and copy to MultiTract/
% use the same name as .json file for the folder
[~,symptomTractsfileName] = fileparts(symptomTractsfile);
ea_connectome_from_pathways(obj, symptomTractsfileName, size(vals,1), size(vals,2), negative_vals)

end
