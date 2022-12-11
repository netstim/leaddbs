function ea_file2json(fname_in,fname_out)

legacy_modalities = {'anat_t1','anat_t2','anat_pd','postop_ct','postop_tra','postop_cor','postop_sag','anat_fgatir','fa','dti','dti','dti','anat_t2star','anat_swi'};
bids_modalities = {'T1w','T2w','PDw','CT','ax','cor','sag','FGATIR','fa','dwi','dwi.bval','dwi.bvec','T2starw','SWI'};
rawdata_containers = containers.Map(legacy_modalities,bids_modalities);
opt.FileName = fname_out;
[filepath,filename,~] = fileparts(fname_in);
json_mat = struct();
%function to convert mat files and text
    if endsWith(fname_in,'.mat')
        %special case for legacy dataset  
        %read the input mat
        %dealing with coregistration
        if strcmp(filename,'ea_coreg_approved')
            [coreg_filepath,~,~] = fileparts(fname_in);
            input_mat = load(fname_in);
            coreg_fieldnames = fieldnames(input_mat);
            for i=1:length(coreg_fieldnames)
                no_of_fieldnames = length(strsplit(coreg_fieldnames{i},'_'));
                if no_of_fieldnames > 2
                    split_fieldname = strsplit(coreg_fieldnames{i},'_');
                    new_fieldname = [split_fieldname{1},'_',split_fieldname{2}];
                    if ismember(new_fieldname,legacy_modalities)
                        modality = rawdata_containers(new_fieldname);
                    elseif contains(new_fieldname, legacy_modalities)
                        %try again because the coreg filename and .json
                        %should be similar
                        match_idx = find(cellfun(@(x) contains(new_fieldname, x), legacy_modalities));
                        if length(match_idx) == 1
                            modality = bids_modalities{match_idx};
                        end
                    else
                        if isvarname(upper(split_fieldname{end}))
                            modality = upper(split_fieldname{end});
                        else
                            modality = upper(split_fieldname{end-1});
                        end
                    end
                    
                else
                    if ismember(coreg_fieldnames{i},legacy_modalities)
                        modality = rawdata_containers(coreg_fieldnames{i});
                   elseif contains(coreg_fieldnames{i}, legacy_modalities)
                        %try again because the coreg filename and .json
                        %should be similar
                        match_idx = find(cellfun(@(x) contains(coreg_fieldnames{i}, x), legacy_modalities));
                        if length(match_idx) == 1
                            modality = bids_modalities{match_idx};
                        end
                    else
                        modality = strsplit(coreg_fieldnames{i},'_');
                        modality = upper(modality{end});
                    end
                    
                end
                flag = 'approval';
                [modality,json_mat] = add_mod(modality,json_mat,flag);
                json_mat.approval.(modality) = input_mat.(coreg_fieldnames{i});
            end
            if exist(fullfile(coreg_filepath,'ea_coregctmethod_applied.mat'),'file')
                temp_mat = load(fullfile(coreg_filepath,'ea_coregctmethod_applied.mat'));
                method_used = generateMethod(temp_mat,'coregct_method_applied');
                modality = 'CT';
                json_mat.method.(modality) = method_used;
            end
            if exist(fullfile(coreg_filepath,'ea_coregmrmethod_applied.mat'),'file')
                temp_mat = load(fullfile(coreg_filepath,'ea_coregmrmethod_applied.mat'));
                method_used = generateMethod(temp_mat,'coregmr_method_applied');
                modality = 'MR';
                json_mat.method.(modality) = method_used;
                temp_fieldname = fieldnames(temp_mat);
                for j=1:length(temp_fieldname)
                    if ~strcmp(temp_fieldname{j},'coregmr_method_applied')
                        no_of_fieldnames = length(strsplit(temp_fieldname{j},'_'));
                        if no_of_fieldnames > 2
                            split_fieldname = strsplit(temp_fieldname{j},'_');
                            new_fieldname = [split_fieldname{1},'_',split_fieldname{2}];
                            if ismember(new_fieldname,legacy_modalities)
                                modality = rawdata_containers(new_fieldname);
                            elseif contains(new_fieldname, legacy_modalities)
                                %try again because the coreg filename and .json
                                %should be similar
                                match_idx = find(cellfun(@(x) contains(new_fieldname, x), legacy_modalities));
                                if length(match_idx) == 1
                                    modality = bids_modalities{match_idx};
                                end
                            else
                                if isvarname(upper(split_fieldname{end}))
                                    modality = upper(split_fieldname{end});
                                else
                                    modality = upper(split_fieldname{end-1});
                                end
                            end
                        else
                            if ismember(temp_fieldname{j},legacy_modalities)
                                modality = rawdata_containers(temp_fieldname{j});
                            elseif contains(temp_fieldname{j}, legacy_modalities)
                                %try again because the coreg filename and .json
                                %should be similar
                                match_idx = find(cellfun(@(x) contains(temp_fieldname{j}, x), legacy_modalities));
                                if length(match_idx) == 1
                                    modality = bids_modalities{match_idx};
                                end
                            else
                                modality = strsplit(temp_fieldname{j},'_');
                                modality = upper(modality{end});
                            end

                        end
                    end
                    flag = 'method';
                    [modality,json_mat] = add_mod(modality,json_mat,flag);
                    json_mat.method.(modality) = temp_mat.(temp_fieldname{j});
                end
            end
            savejson('',json_mat,'method',json_mat.method,opt);
        elseif strcmp(filename,'ea_coregctmethod_applied') && ~exist(fullfile(filepath,'ea_coreg_approved.mat'),'file')
          input_mat = load(fname_in);  
          method_used = generateMethod(input_mat,'coregct_method_applied');
          modality = 'CT';
          json_mat.method.(modality) = method_used;
          savejson('',json_mat,opt);
        elseif strcmp(filename,'ea_coregmrmethod_applied') && ~exist(fullfile(filepath,'ea_coreg_approved.mat'),'file')
          input_mat = load(fname_in);  
          method_used = generateMethod(input_mat,'coregmr_method_applied');
          modality = 'MR';
          json_mat.method.(modality) = method_used;
          savejson('',json_mat,opt);
        %dealing with normalization
        elseif strcmp(filename,'ea_normmethod_applied')
            input_mat = load(fname_in);
            normalize_fieldnames = fieldnames(input_mat);
            for i=1:length(normalize_fieldnames)
                method_used = generateMethod(input_mat,'norm_method_applied');
                json_mat.approval = 1;
                json_mat.method = method_used;
            end
            savejson('',json_mat,opt);
           
        end
    end
 function  method_used = generateMethod(input_mat,modality_field)
     if iscell(modality_field)
         modality_field = modality_field{1};
     end
     if isfield(input_mat,modality_field)
         if ischar(input_mat.(modality_field))
             input_mat.(modality_field) = {input_mat.(modality_field)};
         end
         if strcmp(input_mat.(modality_field){end},'ea_normalize_apply_normalization')
             j = length(input_mat.(modality_field));
             for i=1:length(input_mat.(modality_field))
                 if ~strcmp(input_mat.(modality_field){j},'ea_normalize_apply_normalization')
                     modField = input_mat.(modality_field){j};
                     %need to get the first entry thats not
                     %ea_normalize_apply_normalization
                     break
                 end
                 j = j-1;
             end
         else
             modField = input_mat.(modality_field){end};
         end
         if strcmp(modField,'ANTs') || contains(modField,'_ants')
             method_used = ea_normalize_ants('promt');
         elseif strcmp(modField,'BRAINSFit')
             method_used = 'BRAINSFit (Johnson 2007)';
         elseif strcmp(modField,'FLIRT')
             method_used = 'FLIRT (Jenkinson 2001 & 2002)';
         elseif strcmp(modField,'BBR')
             method_used = 'FLIRT BBR (Greve and Fischl 2009)';
         elseif strcmp(modField,'Hybrid SPM & ANTs')
             method_used = 'Hybrid SPM & ANTs';
         elseif strcmp(modField,'Hybrid SPM & BRAINSFIT')
             method_used = 'Hybrid SPM & BRAINSFIT';
         elseif strcmp(modField,'SPM') || contains(modField,'_spm')
             method_used = 'SPM (Friston 2007)';
         elseif strcmp(modField,'Hybrid SPM & FLIRT')
             method_used = 'Hybrid SPM & FLIRT';
         elseif strcmp(modField,'Schoenecker 2009')
             method_used = 'Three-step affine normalization (ANTs; Schonecker 2009)';
         else
             method_used = '';
         end
     else
         method_used = '';
     end
     return
function [new_mod,json_mat] = add_mod(modality,json_mat,flag)
if ~isempty(fieldnames(json_mat))
    if isfield(json_mat.(flag),modality)
        replaced_old_modality = [modality,num2str(1)];
        %replace old modality name with the new
        json_mat.(flag).(replaced_old_modality) = json_mat.(flag).(modality);
        json_mat.(flag) = rmfield(json_mat.(flag),modality);
        new_mod = [modality,num2str(2)];
    elseif isfield(json_mat.(flag),[modality,num2str(1)])
        suffix = 2;
        modality_to_check = [modality,num2str(suffix)];
        while isfield(json_mat.(flag),modality_to_check)
            suffix = suffix + 1;
            modality_to_check = [modality,num2str(suffix)];
        end
        new_mod = modality_to_check;
    else
        new_mod = modality;
    end
else
    new_mod = modality;
end
return