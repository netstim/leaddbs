function file2json(fname_in,fname_out)

%function to convert mat files and text
    if endsWith(fname_in,'.mat')
        %special case for legacy dataset  
        %read the input mat
        disp(fname_in);        
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
                    if isfield(temp_mat,'coregct_method_applied')
                        if contains(temp_mat.coregct_method_applied{1},'_ants')
                            method_used = ea_normalize_ants('promt');
                        elseif contains(temp_mat.coregct_method_applied{1},'BRAINSFit')
                            method_used = 'BRAINSFit (Johnson 2007)';
                        elseif contains(temp_mat.coregct_method_applied{1},'FLIRT')
                            method_used = 'FLIRT (Jenkinson 2001 & 2002)';
                        elseif contains(temp_mat.coregct_nethod_applied{1},'_spm')
                            method_used = legacyMethods.containers('_spm');
                            json_mat.method.(modality) = method_used;
                        else
                            json_mat.method.(modality) = '';
                        end
                    elseif exist(fullfile(filepath,'ea_coregmrmethod_applied.mat'),'file')
                        temp_mat = load(fullfile(filepath,'ea_coregmrmethod_applied.mat'));
                        if isfield(temp_mat,'coregmr_method_applied')
                            if contains(temp_mat.coregmr_method_applied{1},'_ants')
                                method_used = ea_normalize_ants('promt');
                            elseif contains(temp_mat.coregmr_method_applied{1},'BRAINSFit')
                                method_used = 'BRAINSFit (Johnson 2007)';
                            elseif contains(temp_mat.coregmr_method_applied{1},'FLIRT')
                                method_used = 'FLIRT (Jenkinson 2001 & 2002)';
                            elseif contains(temp_mat.coregmr_method_applied{1},'BBR')
                                method_used = 'FLIRT BBR (Greve and Fischl 2009)';
                            elseif contains(temp_mat.coregmr_method_applied{1},'HybridSPMANTs')
                                method_used = 'Hybrid SPM & ANTs';
                            elseif contains(temp_mat.coregmr_method_applied{1},'HybridSPMBRAINSFIT')
                                method_used = 'Hybrid SPM & BRAINSFIT';
                            elseif contains(temp_mat.coregmr_method_applied{1},'_spm')
                                method_used = 'SPM (Friston 2007)';
                            elseif contains(temp_mat.coregmr_method_applied{1},'HybridSPMFLIRT')
                                method_used = 'Hybrid SPM & FLIRT';
                            end
                            json_mat.method.(modality) = method_used;
                        else
                            json_mat.method.(modality) = '';
                        end
                    end
                end
                encodejson = jsonencode(json_mat);
                json_fid = fopen(fname_out,'w');
            end
        end
        if contains(fname_in,'ea_normmethod_applied.mat')
            input_mat = load(fname_in);
            coreg_fieldnames = fieldnames(input_mat);
            for i=1:length(coreg_fieldnames)
                modality = strsplit(coreg_fieldnames{i},'_');
                modality = upper(modality{end});
                json_mat.approval.(modality) = input_mat.(coreg_fieldnames{i});
            end
        end
        
        %%dealing with normalization
        
    elseif endsWith(fname_in,'.txt')
        fid = fopen(fname_in,'rt');
        json_fid = fopen(fname_out,'w');
        text_cell = textscan(fid, '%s','Delimiter','');
        for text = 1:(length(text_cell{1,1})-1)
            if strcmp(text_cell{1,1}(text),'***')
                    temp_var = strsplit(text_cell{1,1}{text+1},':');
                    method_var = strtrim(temp_var{end});
                    S.method.(method_var) = text_cell{1,1}{text+3};
            end
        end
        if ~exist('S','var')
            try
                %try a type of manual resolution - not sure if it belongs to the catch statement?
                temp_var = strsplit(text_cell{1,1}{1},':');
                method_var = strtrim(temp_var{end});
                S.method.(method_var) = text_cell{1,1}{3};
            catch
                disp(['You might have to convert this file to .json manually: ' fname_in]);
            end
        end
        if exist('S','var')
            encodejson = jsonencode(S);
        end
    end
    if exist('encodejson','var') && exist('json_fid','var')
        fprintf(json_fid,encodejson);
    end
    end