function ea_file2json(fname_in,fname_out)

legacy_modalities = {'anat_t1','anat_t2','anat_pd','postop_ct','postop_tra','postop_cor','postop_sag','anat_fgatir','fa','dti','dti','dti','anat_t2star'};
bids_modalities = {'T1w','T2w','PDw','CT','ax','cor','sag','FGATIR','fa','dwi','dwi.bval','dwi.bvec','T2starw'};
rawdata_containers = containers.Map(legacy_modalities,bids_modalities);
opt.FileName = fname_out;
[filepath,filename,~] = fileparts(fname_in);
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
                if ismember(coreg_fieldnames{i},legacy_modalities)
                    modality = rawdata_containers(coreg_fieldnames{i});
                else
                    modality = strsplit(coreg_fieldnames{i},'_');
                    modality = upper(modality{end});
                end
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
                        if ismember(temp_fieldname{j},legacy_modalities)
                            modality = rawdata_containers(temp_fieldname{j});
                        else
                            modality = strsplit(temp_fieldname{j},'_');
                            modality = upper(modality{end});
                        end
                        json_mat.method.(modality) = temp_mat.(temp_fieldname{j});
                    end
                end
            end
            savejson('',json_mat,'method',json_mat.method,opt);
        elseif strcmp(filename,'ea_coregctmethod_applied') && ~exist(fullfile(filepath,'ea_coreg_approved.mat'),'file')
          input_mat = load(fname_in);  
          method_used = generateMethod(input_mat,'coregct_method_applied');
          modality = 'CT';
          json_mat.method.(modality) = method_used;
          savejson('',json_mat,'method',json_mat.method,opt);
        elseif strcmp(filename,'ea_coregmrmethod_applied') && ~exist(fullfile(filepath,'ea_coreg_approved.mat'),'file')
          input_mat = load(fname_in);  
          method_used = generateMethod(input_mat,'coregmr_method_applied');
          modality = 'MR';
          json_mat.method.(modality) = method_used;
          savejson('',json_mat,'method',json_mat.method,opt);
        %dealing with normalization
        elseif strcmp(filename,'ea_normmethod_applied')
            input_mat = load(fname_in);
            normalize_fieldnames = fieldnames(input_mat);
            for i=1:length(normalize_fieldnames)
                method_used = generateMethod(input_mat,'norm_method_applied');
                json_mat.approval = 1;
                json_mat.method = method_used;
            end
            savejson('',json_mat,'method',json_mat.method,'approval',json_mat.approval,opt);
           
        end
    end
 function  method_used = generateMethod(input_mat,modality_field)
     if iscell(modality_field)
         modality_field = modality_field{1};
     end
     if isfield(input_mat,modality_field) && iscell(modality_field)
         if strcmp(input_mat.(modality_field){end},'ANTs') || contains(input_mat.(modality_field){end},'_ants')
             method_used = ea_normalize_ants('promt');
         elseif strcmp(input_mat.(modality_field){end},'BRAINSFit')
             method_used = 'BRAINSFit (Johnson 2007)';
         elseif strcmp(input_mat.(modality_field){end},'FLIRT')
             method_used = 'FLIRT (Jenkinson 2001 & 2002)';
         elseif strcmp(input_mat.(modality_field){end},'BBR')
             method_used = 'FLIRT BBR (Greve and Fischl 2009)';
         elseif strcmp(input_mat.(modality_field){end},'Hybrid SPM & ANTs')
             method_used = 'Hybrid SPM & ANTs';
         elseif strcmp(input_mat.(modality_field){end},'Hybrid SPM & BRAINSFIT')
             method_used = 'Hybrid SPM & BRAINSFIT';
         elseif strcmp(input_mat.(modality_field){end},'SPM') || contains(input_mat.(modality_field){end},'_spm')
             method_used = 'SPM (Friston 2007)';
         elseif strcmp(input_mat.(modality_field){end},'Hybrid SPM & FLIRT')
             method_used = 'Hybrid SPM & FLIRT';
         else
             method_used = '';
         end
     else
         method_used = '';
     end
     return
