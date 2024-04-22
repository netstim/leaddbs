function ea_get_native_Eproj(tractset,stim_ID)
% Comptute E-field metrics on fibers warped into native space.
% The results are stored as
% patient_folder/miscellaneous/connnectomeName/stim_ID/E_peak.mat (values for all fibers in the connectome!)
% By Butenko and Roediger, konstantinmgtu@gmail.com

arguments
    tractset       % fiber filtering object
    stim_ID        % when using PseudoM, provide stimulation folder full(!) name
end

%obj.connectome = 'PPU_rh_downsampled_by_4';

% use connectome with all pathways combined
if tractset.multi_pathways == 1
    % check if merged_pathways is in fibfiltering folder
    [filepath,~,~] = fileparts(tractset.analysispath);
    merged_connectome = [filepath,filesep,tractset.connectome,filesep,'merged_pathways.mat'];
    if ~isfile(merged_connectome)
        % else check if it is in the original lead-group folder
        [filepath,~,~] = fileparts(tractset.leadgroup);
        merged_connectome = [filepath,filesep,tractset.connectome,filesep,'merged_pathways.mat'];
        if ~isfile(merged_connectome)
            % or if it is in another lead-group folder (where fibfiltering file is)
            [filepath,~,~] = fileparts(tractset.analysispath);
            [filepath,~,~] = fileparts(filepath);
            merged_connectome = [filepath,filesep,tractset.connectome,filesep,'merged_pathways.mat'];
        end
    end
else
    merged_connectome = [ea_getconnectomebase('dMRI'), tractset.connectome, filesep, 'data.mat'];
end


% find where VATs are
if ~isfield(tractset.M,'pseudoM')
    stim_space = ['/native/gs_',tractset.M.guid]; 
else
    % custom case
    stim_space = ['/native/',stim_ID];
end

for pt_i = 1:size(tractset.M.patient.list,1)

    stim_folder = strcat(tractset.M.patient.list{pt_i},filesep,'stimulations',stim_space);

    for side = 1:2
    
        switch side
            case 1
                side_suffix = '_rh';
            case 2
                side_suffix = '_lh';
        end

        % OSS-DBS format
        result_folder = strcat(stim_folder,'/Results',side_suffix);
       
        if isfolder(result_folder)
            % get the solution(s)
            mycsv = dir(fullfile(result_folder,'E_field_Lattice*'));
    
            % create 4D nii (4-th dimension is for E-field components and magnitude)
            for field_i = 1:length(mycsv)
                ea_get_4Dfield_from_csv(mycsv(field_i).folder, mycsv(field_i).name, result_folder)
            end
    
            myFields = dir(fullfile(result_folder,'/4D_E_field_Lattice*')); % gets all mat files in struct
            % compute projection of the E-fields onto the fibers
            for field_i = 1:length(myFields)
                e_field_file = fullfile(myFields(field_i).folder, myFields(field_i).name);
                ea_get_E_field_along_fibers(tractset.M.patient.list{pt_i},['gs_',tractset.M.guid],e_field_file, merged_connectome, side_suffix,tractset.calcthreshold)
            end
        else
            [~,pt_label,~] = fileparts(tractset.M.patient.list{pt_i});
            fprintf("Missing stimulation for %s, %s side \n",pt_label,side_suffix)
        end

        disp("_______________________________")
    end

end

