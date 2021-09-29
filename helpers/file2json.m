function file2json(fname_in,fname_out)

%function to convert mat files and text
    if endsWith(fname_in,'.mat')
        %special case for legacy dataset  
        %read the input mat
        %dealing with coregistration
        if contains(fname_in,'coreg_approved.mat')
            [filepath,filename,ext] = fileparts(fname_in);
            input_mat = load(fname_in);
            coreg_fieldnames = fieldnames(input_mat);
            for i=1:length(coreg_fieldnames)
                modality = strsplit(coreg_fieldnames{i},'_');
                modality = upper(modality{end});
                json_mat.approval.(modality) = input_mat.(coreg_fieldnames{i});
                if exist(fullfile(filepath,'ea_coregctmethod_applied.mat'),'file')
                    temp_mat = load(fullfile(filepath,'ea_coregctmethod_applied.mat'));
                    method_used = generateMethod(temp_mat,'coregct_method_applied');
                elseif exist(fullfile(filepath,'ea_coregmrmethod_applied.mat'),'file')
                    temp_mat = load(fullfile(filepath,'ea_coregmrmethod_applied.mat'));
                    method_used = generateMethod(temp_mat,'coregmr_method_applied');
                end
                json_mat.method.(modality) = method_used;
            end
            encodejson = jsonencode(json_mat);
            json_fid = fopen(fname_out,'w');
            fprintf(json_fid,encodejson);
        
        %dealing with normalization
        elseif contains(fname_in,'ea_normmethod_applied.mat')
            input_mat = load(fname_in);
            normalize_fieldnames = fieldnames(input_mat);
            for i=1:length(normalize_fieldnames)
                method_used = generateMethod(input_mat,'norm_method_applied');
                json_mat.approval = 1;
                json_mat.method = method_used;
            end
            encodejson = jsonencode(json_mat);
            json_fid = fopen(fname_out,'w');
            fprintf(json_fid,encodejson);
        end
    end
 function  method_used = generateMethod(input_mat,modality_field)
     if iscell(modality_field)
         modality_field = modality_field{1};
     end
     if isfield(input_mat,modality_field)
         if contains(input_mat.(modality_field),'_ants')
             method_used = ea_normalize_ants('promt');
         elseif contains(input_mat.(modality_field),'BRAINSFit')
             method_used = 'BRAINSFit (Johnson 2007)';
         elseif contains(input_mat.(modality_field),'FLIRT')
             method_used = 'FLIRT (Jenkinson 2001 & 2002)';
         elseif contains(input_mat.(modality_field),'BBR')
             method_used = 'FLIRT BBR (Greve and Fischl 2009)';
         elseif contains(input_mat.(modality_field),'HybridSPMANTs')
             method_used = 'Hybrid SPM & ANTs';
         elseif contains(input_mat.(modality_field),'HybridSPMBRAINSFIT')
             method_used = 'Hybrid SPM & BRAINSFIT';
         elseif contains(input_mat.(modality_field),'_spm')
             method_used = 'SPM (Friston 2007)';
         elseif contains(input_mat.(modality_field),'HybridSPMFLIRT')
             method_used = 'Hybrid SPM & FLIRT';
         end
     else
         method_used = '';
     end
     return
%%dealing with normalization
        
%     elseif endsWith(fname_in,'.txt')
%         fid = fopen(fname_in,'rt');
%         json_fid = fopen(fname_out,'w');
%         text_cell = textscan(fid, '%s','Delimiter','');
%         for text = 1:(length(text_cell{1,1})-1)
%             if strcmp(text_cell{1,1}(text),'***')
%                     temp_var = strsplit(text_cell{1,1}{text+1},':');
%                     method_var = strtrim(temp_var{end});
%                     S.method.(method_var) = text_cell{1,1}{text+3};
%             end
%         end
%         if ~exist('S','var')
%             try
%                 %try a type of manual resolution - not sure if it belongs to the catch statement?
%                 temp_var = strsplit(text_cell{1,1}{1},':');
%                 method_var = strtrim(temp_var{end});
%                 S.method.(method_var) = text_cell{1,1}{3};
%             catch
%                 disp(['You might have to convert this file to .json manually: ' fname_in]);
%             end
%         end
%         if exist('S','var')
%             encodejson = jsonencode(S);
%         end
%     end
%     if exist('encodejson','var') && exist('json_fid','var')
%         fprintf(json_fid,encodejson);
%     end
%     end