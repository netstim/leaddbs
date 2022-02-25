function results = ea_regenerateRawjson(source_dir)
 if ~iscell(source_dir)
     source_dir ={source_dir};
 end
 rawdata_folder = fullfile(source_dir,'rawdata');
 for i=1:length(rawdata_folder)
     subj_name = dir_without_dots(rawdata_folder{i});
     subj_name = subj_name([subj_name.isdir]);
     subj_name = subj_name.name;
     prefs_dir = fullfile(source_dir{i},'derivatives','leaddbs',subj_name,'prefs');
     if ~exist(prefs_dir,'dir')
         disp("preparing prefs directory...")
         mkdir(prefs_dir)
     end
     raw_preop_dir = fullfile(rawdata_folder{i},subj_name,'ses-preop','anat');
     raw_postop_dir = fullfile(rawdata_folder{i},subj_name,'ses-postop','anat');
     opt.FileName = fullfile(prefs_dir,[subj_name,'_','desc-rawimages.json']);
     
     preop_files = dir_without_dots(raw_preop_dir);
     preop_files = {preop_files.name};
     
     if ~isempty(preop_files)
         for j=1:length(preop_files)
             json_val = regexprep(preop_files{j},'(.nii)|(.gz)','');
             if contains(preop_files{j},'acq-')
                 temp_tag = strsplit(json_val,'-');
                 modality_str = temp_tag{end};
                 rawdata_fieldname = modality_str;
                 anat_files_selected.preop.anat.(rawdata_fieldname) = json_val;
             else
                 temp_tag = strsplit(json_val,'preop_');
                 rawdata_fieldname = temp_tag{end};
                 anat_files_selected.preop.anat.(rawdata_fieldname) = json_val;
             end
         end
     end
     postop_files = dir_without_dots(raw_postop_dir);
     postop_files = {postop_files.name};
     
     for k=1:length(postop_files)
         json_val = regexprep(postop_files{k},'(.nii)|(.gz)','');
         if contains(postop_files{k},'ct','IgnoreCase',true)
             rawdata_fieldname = 'CT';
             anat_files_selected.postop.anat.(rawdata_fieldname) = json_val;
         else
             temp_tag = strsplit(json_val,'-');
             rawdata_fieldname = temp_tag{end};
             anat_files_selected.postop.anat.(rawdata_fieldname) = json_val;
             
         end
     end
     disp(['Generated rawimages.json. Saving as:',opt.FileName]);
     savejson('',anat_files_selected,opt);
     
 end
end