function ea_write_scores(M,nummat,matname,new_var,val_to_rm)
%transfers to the group folder

if isfield(M,'clinical')
    for pt=1:length(M.patient.list)
        score_dir = [M.patient.list{pt},'/clinical','/clinical_scores.mat'];
        if exist(fullfile(M.patient.list{pt},'clinical','clinical_scores.mat'),'file')
            load(fullfile(M.patient.list{pt},'clinical','clinical_scores.mat'));
            if strcmp(fieldnames(scores),'Motor_MDSUPDRS')
                score_type = 'Motor_MDSUPDRS';
            elseif strcmp(fieldnames(scores),'Motor_UPDRS')
                score_type = 'Motor_UPDRS';
            elseif strcmp(fieldnames(scores),'BDI')
                score_type = 'BDI';
            else
                score_cell = fieldnames(scores);
                score_type = score_cell{1};
            end
            if new_var == 0 % changing the old variable
                matcell = strsplit(matname,'-');
                score_type = matcell{1};
                postop_flag = matcell{2};
                fieldname = matcell{3};
                val_name = matcell{4};
                scores.(score_type).(postop_flag).(fieldname).(val_name) = nummat(pt);
                save(score_dir,'scores')
                if pt==1
                    disp("Storing new clinical scores into patient folders")
                end
            elseif new_var == 1 %ading a new variable
                  postopid = fieldnames(scores.(score_type));
                  postop_flag = postopid{1,1};
                  scores.(score_type).(postop_flag).(matname).value = nummat(pt);
                  save(score_dir,'scores')
                  if pt==1
                        disp("Storing modified clinical stores into patient folders")
                  end  
            end
            if ~strcmp(val_to_rm,'')            
                str_to_rmv = strsplit(val_to_rm{1},'-');
                postop_flag = str_to_rmv{1};
                if isfield(scores.(score_type).(postop_flag),(str_to_rmv{2}))
                    field = str_to_rmv{2};
                    scores.(score_type).(postop_flag) = rmfield(scores.(score_type).(postop_flag),field);
                    save(score_dir,'scores')
                    if pt==1
                        disp("Removing clinical scores from patient directory")
                    end
                end
            end
        end
    end
end