function ea_write_scores(M,nummat,matname,new_var,val_to_rm)
%transfers to the group folder

if isfield(M,'clinical')
    for pt=1:length(M.patient.list)
        score_dir = [M.patient.list{pt},'/clinical','/clinical_scores.mat'];
        if exist(fullfile(M.patient.list{pt},'clinical','clinical_scores.mat'),'file')
            load(fullfile(M.patient.list{pt},'clinical','clinical_scores.mat'));
            score_cell = fieldnames(scores);
            
        end
        if new_var == 0 % changing the old variable
            matcell = strsplit(matname,'-');
            score_type = matcell{1};
            if strcmp(score_type,'Motor_Mixed')
                score_type = 'Motor_UPDRS';
                if ~isfield(scores,score_type)
                    score_type = 'Motor_MDSUPDRS';
                end
            end
            postop_flag = matcell{2};
            fieldname = matcell{3};
            val_name = matcell{4};
            scores.(score_type).(postop_flag).(fieldname).(val_name) = nummat(pt);
            save(score_dir,'scores')
            if pt==1
                disp("Storing modified clinical scores into patient folders")
            end
        elseif new_var == 1 %ading a new variable
            score_type = 'NEW';
            postop_flag = 'Custom';
            if contains(matname,' ')
                matname = strrep(matname,' ','');
            end
            scores.(score_type).(postop_flag).(matname).value = nummat(pt);
            save(score_dir,'scores')
            if pt==1
                disp("Storing new clinical stores into patient folders")
            end
        end
        if ~strcmp(val_to_rm,'')
            str_to_rmv = strsplit(val_to_rm{1},'-');
            score_type = str_to_rmv{1};
            postop_flag = str_to_rmv{2};
            sub_score = str_to_rmv{3};
            if strcmp(score_type,'Motor_Mixed')
                score_type = 'Motor_UPDRS';
                if ~isfield(scores,score_type)
                    score_type = 'Motor_MDSUPDRS';
                end
            end
            if isfield(scores.(score_type).(postop_flag),(sub_score))
                field_type = str_to_rmv{4};
                scores.(score_type).(postop_flag).(sub_score) = rmfield(scores.(score_type).(postop_flag).(sub_score),field_type);
                save(score_dir,'scores')
                if pt==1
                    disp("Removing clinical scores from patient directory")
                end
            
                
            end
        end
    end
end
end