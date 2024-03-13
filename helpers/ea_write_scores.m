function ea_write_scores(M,nummat,matname,new_var,val_to_rm)
%transfers to the group folder

if isfield(M,'clinical')
    for pt=1:length(M.patient.list)
        [~,subj_id,~] = fileparts(M.patient.list{pt});
        guid = ['gs_' M.guid];
        score_file = fullfile(M.patient.list{pt},'clinical',guid,[subj_id,'_desc-clinicalScores.mat']);
        if exist(score_file,'file')
            load(score_file);
        end
        if new_var == 0 % changing the old variable
            matcell = strsplit(matname,'-');
            score_type = matcell{1};
            postop_flag = matcell{2};
            val_name = matcell{3};
            if pt==1
                disp("Storing modified clinical scores into patient folders")
            end
        elseif new_var == 1 %ading a new variable
            score_type = 'newvar';
            postop_flag = 'postop';
            val_name = 'improvments';
            if contains(matname,' ')
                matname = strrep(matname,' ','');
            end
            if pt==1
                disp("Storing new clinical stores into patient folders")
            end
        end
        clinical.(guid).scores.(postop_flag).(score_type).(val_name) = nummat(pt);
        if ~strcmp(val_to_rm,'')
            matcell = strsplit(val_to_rm{1},'-');
            score_type = matcell{1};
            postop_flag = matcell{2};
            val_name = matcell{3};
            if isfield(clinical.(guid).scores.(postop_flag),score_type)
                clinical.(guid).scores.(postop_flag).(score_type) = rmfield(clinical.(guid).scores.(postop_flag).(score_type),val_name);
                if pt==1
                    disp("Removing clinical scores from patient directory")
                end
            end
        end
        save(score_file,'clinical')
    end
end
end