function [mni_files,native_files] = vta_walkpath(new_path,pipeline)


switch pipeline
    case 'stimulations'
        mni_dir = fullfile(new_path,pipeline,'MNI_ICBM_2009b_NLIN_ASYM');
        native_dir = fullfile(new_path,pipeline,'native');
        if exist(mni_dir,'dir')
            this_folder = dir_without_dots(mni_dir);
            this_folder_names = {this_folder.name};
            file_indx = 1;
            for folder_names = 1:length(this_folder_names)
                files_in_this_folder = dir_without_dots(fullfile(mni_dir,this_folder_names{folder_names}));
                mni_files{file_indx} = {files_in_this_folder.name};
                for mni_file = 1:length(mni_files{1,file_indx})
                    mni_files{1,file_indx}{1,mni_file} = fullfile(mni_dir,this_folder_names{folder_names},mni_files{1,file_indx}{1,mni_file});
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
                for native_file = 1:length(native_files{1,file_indx})
                    native_files{1,file_indx}{1,native_file} = fullfile(native_dir,this_folder_names{folder_names},native_files{1,file_indx}{1,native_file});
                end
                file_indx = file_indx+1;
                %all_files contains all the names
                %of the vta files present. Now you
                %can rename them.
            end
        end
        
    case 'current_headmodel'
      mni_dir = fullfile(new_path,pipeline,'MNI_ICBM_2009b_NLIN_ASYM');
      native_dir = fullfile(new_path,pipeline,'native');  
      if exist(mni_dir,'dir')
          this_folder = dir_without_dots(mni_dir);
          mni_files = {this_folder.name};
          for i=1:length(mni_files)
              mni_files{i} = fullfile(mni_dir,mni_files{i});
          end
      end
      if exist(native_dir,'dir')
          this_folder = dir_without_dots(native_dir);
          native_files = {this_folder.name};
          for i=1:length(native_files)
              native_files{i} = fullfile(native_dir,native_files{i});
          end
      end
end
return