function ea_file2json(fname_in,fname_out)

legacy_modalities = {'anat_t1','anat_t2','anat_pd','postop_ct','postop_tra','postop_cor','postop_sag','anat_fgatir','anat_flair','fa','dti','dti','dti','anat_t2star','anat_swi'};
bids_modalities = {'T1w','T2w','PDw','CT','ax_MRI','cor_MRI','sag_MRI','FGATIR','FLAIR','fa','dwi','dwi.bval','dwi.bvec','T2starw','SWI'};
rawdata_containers = containers.Map(bids_modalities,legacy_modalities);
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
            [coregDir,~,~] = fileparts(fname_out);
            coregDir = strrep(coregDir,'log','anat');
            input_mat = load(fname_in);
            coreg_fieldnames = fieldnames(input_mat);
            %determine files inside the coreg files
            coregdir_filenames = ea_regexpdir(coregDir,'.*.nii',0,'f');
            for i=1:length(coregdir_filenames)
                if contains(coregdir_filenames{i},'acq-')
                   
                    tmpmod = strsplit(coregdir_filenames{i},'acq-');
                    mod = strrep(tmpmod{end},'.nii','');
                    tmpmod = regexprep(mod, '(_)?(iso|sag|ax|cor)\d?(_)?', '');
                    try
                        %now find this in the coregapproved file
                        if contains(mod,'_MRI')
                            legacymod = rawdata_containers(mod);
                        else
                            legacymod = rawdata_containers(tmpmod);
                        end
                    catch
                        testind = cellfun(@(x)isequal(x,tmpmod),bids_modalities);
                        if any(testind)
                            legacymod = legacy_modalities{testind};
                        else
                            legacymod = lower(tmpmod);
                        end
                    end
                    idx = cellfun(@(x)contains(x,legacymod),coreg_fieldnames);
                    if ~isempty(find(idx,1))
                        if length(find(idx)) > 1
                            acq = strsplit(mod,'_');
                            mod_woacq = acq{1};
                            digit = regexp(mod_woacq, '\d', 'match');
                            if ~isempty(digit)
                                digit = digit{1};
                                corridx = find(idx);
                                idx = corridx(str2double(digit));
                                json_mat.approval.(mod) = input_mat.(coreg_fieldnames{idx});
                            else
                                idx = cellfun(@(x)strcmp(x,legacymod),coreg_fieldnames);
                                if any(idx)
                                %quite difficult to find extreme cases. for
                                %now, make their approval 0
                                    json_mat.approval.(mod) = input_mat.(coreg_fieldnames{idx});
                                else
                                     json_mat.approval.(mod) = 0;
                                     warning("One or more elements of the coreg fieldname were not transformed. Please review your desc_coregmethod.json file inside the derivatives/leaddbs/coregistration/log directory...")
                                end
                                
                            end                        
                        else
                            json_mat.approval.(mod) = input_mat.(coreg_fieldnames{idx});
                        end
                    else
                        json_mat.approval.(mod) = 0;
                        warning("One or more elements of the coreg fieldname were not transformed. Please review your desc_coregmethod.json file inside the derivatives/leaddbs/coregistration/log directory...")
                    end
                end
            end
            savejson('',json_mat,'approval',json_mat.approval,opt);   
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
 function method_used = generateMethod(input_mat,modality_field)
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
             warning("We could not identify the normalization method used. Please edit it manually.\n" + ...
                 " You will find the file under derivatives/leaddbs/normalization/log!");
         end
     else
         method_used = '';
         warning("We could not identify the normalization method used. Please edit it manually.\n" + ...
                 " You will find the file under derivatives/leaddbs/normalization/log!");
     end
     return
