function [mni_files,native_files,derivatives_cell,mni_model_names,native_model_names] = ea_vta_walkpath(source_patient,new_path,pipeline,derivatives_cell)
[~,~,~,~,~,~,~,~,~,~,~,lead_mapper] = ea_create_bids_mapping();
%subjID = regexp(new_path,'[sub-].*[Aa-Zz][0-9]$','match');
extract_subjID_str = strsplit(new_path,'leaddbs/');
subjID = extract_subjID_str{2};
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
            old_mni_dir = fullfile(source_patient,pipeline,'MNI_ICBM_2009b_NLIN_ASYM');
            this_folder = dir_without_dots(mni_dir);
            this_folder_names = {this_folder.name};
            file_indx = 1;
            for folder_names = 1:length(this_folder_names)
                files_in_this_folder = dir_without_dots(fullfile(mni_dir,this_folder_names{folder_names}));
                % files_in_this_folder = files_in_this_folder(~files_in_this_folder.isdir);
                mni_files{file_indx} = {files_in_this_folder.name};
                mni_model_name = add_model(fullfile(mni_dir,this_folder_names{folder_names}));
                mni_model_names{folder_names} = mni_model_name;
                for mni_file = 1:length(mni_files{1,file_indx})
                    if isfolder(fullfile(mni_dir,this_folder_names{folder_names},mni_files{1,file_indx}{1,mni_file}))
                        possible_connectome_folder = fullfile(mni_dir,this_folder_names{folder_names},mni_files{1,file_indx}{1,mni_file});
                        possible_MRI_file = ea_regexpdir(possible_connectome_folder,'.*vat_seed_compound_[df]MRI.*',0);
                        if ~isempty(possible_MRI_file)
                            bids_connectome_name = ea_getConnLabel(mni_files{1,file_indx}{1,mni_file});
                            stimulation_folder = fullfile(mni_dir,this_folder_names{folder_names});
                            connectome_folder = fullfile(mni_dir,this_folder_names{folder_names},mni_files{1,file_indx}{1,mni_file});
                            %create a struct of all the properties
                            %associated with this file
                            tag_struct.subjID = subjID;
                            tag_struct.modeltag = mni_model_name;
                            tag_struct.conntag = bids_connectome_name;
                            generate_bidsConnectome_name(stimulation_folder,connectome_folder,lead_mapper,tag_struct)
                        end
                    else
                        derivatives_cell{end+1,1} = fullfile(old_mni_dir,mni_files{1,file_indx}{1,mni_file});
                        mni_files{1,file_indx}{1,mni_file} = fullfile(mni_dir,this_folder_names{folder_names},mni_files{1,file_indx}{1,mni_file});
                        derivatives_cell{end,2} = mni_files{1,file_indx}{1,mni_file};
                    end
                end
                file_indx = file_indx+1;
                %all_files contains all the names
                %of the vta files present. Now you
                %can rename them.
            end

        end
        if exist('native_dir','dir')
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
      
end
return
end

function model_name = add_model(stimFolder)
  stimParams = ea_regexpdir(stimFolder, 'stimparameters\.mat$', 0);
  if ~isempty(stimParams)
    load(stimParams{1},'S')
    model_name = ea_simModel2Label(S.model);
  else
      ea_cprintf('CmdWinWarnings', 'Missing stimparameters under %s\nSet to SimBio model by default, please check manually.\n', stimFolder);
      model_name = 'simbio';
  end
  return
end

function generate_bidsConnectome_name(mni_folder,connectome_folder,lead_mapper,tag_struct)
   %connectome_folder should be the full path of the connectome files,
   %mni_folder =
   %fullfile(new_path,pipeline,'MNI152NLinAsym','stimulations',stimulation
   %name)
   mapper_output_files = dir_without_dots(connectome_folder);
   mapper_output_files = {mapper_output_files(~[mapper_output_files.isdir]).name};
   
   for mapper_file = 1:length(mapper_output_files)
       matching_file = regexprep(mapper_output_files{mapper_file},'(fl_|efield_|efield_gauss)', '');
       matching_file = ['_' matching_file];
       if ismember(matching_file,lead_mapper{:,1})
           if contains(mapper_output_files{mapper_file},'efield')
               tag_struct.simtag = 'efield';
           elseif contains(mapper_output_files{mapper_file},'efield_gauss')
               tag_struct.simtag = 'efieldgauss';
           else %doesn't contain either
               tag_struct.simtag = 'binary';
           end
           extract_old_hemisdesc_str = strsplit(mapper_output_files{mapper_file},'_');
           
           indx = cellfun(@(x)strcmp(x,matching_file),lead_mapper{:,1});
           bids_name = lead_mapper{1,2}{indx};
           %replace hemidesc tag
           if strcmp(extract_old_hemisdesc_str{1},'fl')
               bids_name = strrep(bids_name,'flippedTag','flipped');
           else
               bids_name = strrep(bids_name,'_hemidesc-flippedTag','');
           end
           %replace sim tag
           bids_name = strrep(bids_name,'simTag',tag_struct.simtag);
           %replace model tag
           bids_name = strrep(bids_name,'modelTag',tag_struct.modeltag);
           %replace connectome tag
           bids_name = strrep(bids_name,'conTag',tag_struct.conntag);
           copyfile(fullfile(connectome_folder,mapper_output_files{mapper_file}),fullfile(mni_folder,mapper_output_files{mapper_file}));
           movefile(fullfile(mni_folder,mapper_output_files{mapper_file}),fullfile(mni_folder,[tag_struct.subjID,'_',bids_name]));
           fprintf('Renaming file %s as %s.\n',mapper_output_files{mapper_file},[tag_struct.subjID,'_',bids_name]);
           
       else
           evalin('base','WARNINGSILENT=1;');
           ea_warning(sprintf('BIDS tag could not be assigned for %s. Please rename manually',mapper_output_files{mapper_file}));
       end
       
   end
   evalin('base','WARNINGSILENT=1;');
   ea_warning(sprintf('Deleting old copy of connectome folder %s. You can find it in the source patient folder if you need.',connectome_folder));
   %ea_delete(fullfile(connectome_folder));
end




