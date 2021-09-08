function [mni_files,native_files] = vta_walkpath(new_path,pipeline)



mni_dir = fullfile(new_path,pipeline,'MNI_ICBM_2009b_NLIN_ASYM');
native_dir = fullfile(new_path,pipeline,'native');
if exist(mni_dir,'dir')
    this_folder = dir_without_dots(mni_dir);
    this_folder_names = {this_folder.name};
    for folder_names = 1:length(this_folder_names)
        files_in_this_folder = dir_without_dots(fullfile(mni_dir,this_folder_names{folder_names}));
        mni_files = {files_in_this_folder.name};
        for mni_file = 1:length(mni_files)
            mni_files{mni_files} = fullfile(mni_dir,this_folder_names{folder_names},mni_files{mni_file});
        end
        %all_files contains all the names
        %of the vta files present. Now you
        %can rename them.
    end
end
if exist(native_dir,'dir')
    this_folder = dir_without_dots(native_dir);
    this_folder_names = {this_folder.name};
    for folder_names = 1:length(this_folder_names)
        files_in_this_folder = dir_without_dots(fullfile(mni_dir,this_folder_names{folder_names}));
        mni_files = {files_in_this_folder.name};
        for mni_file = 1:length(mni_files)
            mni_files{mni_files} = fullfile(mni_dir,this_folder_names{folder_names},mni_files{mni_file});
        end
        %all_files contains all the names
        %of the vta files present. Now you
        %can rename them.
    end
end