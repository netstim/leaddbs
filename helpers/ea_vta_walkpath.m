function [mni_files,native_files,derivatives_cell,mni_model_names,native_model_names] = ea_vta_walkpath(source_patient,new_path,pipeline,derivatives_cell)

mni_files = {};
native_files = {};
mni_model_names = {};
native_model_names = {};
switch pipeline
    case 'stimulations'
        if ~exist(fullfile(new_path,pipeline),'dir')
           return 
        end
        if exist(fullfile(new_path,pipeline,'MNI152NLin2009bAsym'),'dir')
            mni_dir = fullfile(new_path,pipeline,'MNI152NLin2009bAsym');
            native_dir = fullfile(new_path,pipeline,'native');
            this_folder = dir_without_dots(mni_dir);
            this_folder_names = {this_folder.name};
            file_indx = 1;
            for folder_names = 1:length(this_folder_names)
                files_in_this_folder = dir_without_dots(fullfile(mni_dir,this_folder_names{folder_names}));
                mni_files{file_indx} = {files_in_this_folder.name};
                mni_model_name = add_model(fullfile(source_patient,pipeline,'MNI_ICBM_2009b_NLIN_ASYM',this_folder_names{folder_names}));
                mni_model_names{folder_names} = mni_model_name; 
                for mni_file = 1:length(mni_files{1,file_indx})
                    derivatives_cell{end+1,1} = fullfile(source_patient,pipeline,'MNI_ICBM_2009b_NLIN_ASYM',mni_files{1,file_indx}{1,mni_file});
                    mni_files{1,file_indx}{1,mni_file} = fullfile(mni_dir,this_folder_names{folder_names},mni_files{1,file_indx}{1,mni_file});
                    derivatives_cell{end,2} = mni_files{1,file_indx}{1,mni_file};
                end
                file_indx = file_indx+1;
                %all_files contains all the names
                %of the vta files present. Now you
                %can rename them.
            end            
            
        end
        if exist(native_dir,'dir')
            this_folder = dir_without_dots(native_dir);
            this_folder_names = {this_folder.name};
            file_indx = 1;
            for folder_names = 1:length(this_folder_names)
                files_in_this_folder = dir_without_dots(fullfile(native_dir,this_folder_names{folder_names}));
                native_files{file_indx} = {files_in_this_folder.name};
                native_model_name = add_model(fullfile(source_patient,pipeline,'native',this_folder_names{folder_names}));
                native_model_names{folder_names} = native_model_name; 
                for native_file = 1:length(native_files{1,file_indx})
                    derivatives_cell{end+1,1} = fullfile(source_patient,pipeline,'native',native_files{1,file_indx}{1,native_file});
                    native_files{1,file_indx}{1,native_file} = fullfile(native_dir,this_folder_names{folder_names},native_files{1,file_indx}{1,native_file});
                    derivatives_cell{end,2} = native_files{1,file_indx}{1,native_file};
                end
                file_indx = file_indx+1;
                %all_files contains all the names
                %of the vta files present. Now you
                %can rename them.
            end
        end
  
      if exist(native_dir,'dir')
          this_folder = dir_without_dots(native_dir);
          native_files = {this_folder.name};
          for i=1:length(native_files)
              derivatives_cell{end+1,1} = fullfile(source_patient,pipeline,'native',native_files{i});
              native_files{i} = fullfile(native_dir,native_files{i});
              derivatives_cell{end,2} = native_files{i};
          end
      end
      
end
return
end

function model_name = add_model(stimFolder)
  stimParams = ea_regexpdir(stimFolder, 'stimparameters\.mat$', 0);
  if ~isempty(stimParams)
    load(stimParams{1},'S')
    model_name = ea_simModel2Label(S.model);
  else
      warndlg("Could not detect the model in your parameters file. Please check manually. Adding a default name for now...");
      model_name = 'default';
  end
  return
end

