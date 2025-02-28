function ea_get_Eproj(tractset,vatlist)
% Comptute E-field metrics directly on fibers (optionally in native space).
% The results are stored as
% patient_folder/connectomes/connnectomeName/stim_ID/E_metrics.mat (values for all fibers in the connectome!)
% By Butenko and Roediger, konstantinmgtu@gmail.com

arguments
    tractset       % fiber filtering object
    vatlist        % cell object {pt_N,2} with paths to 4D niftis
end

if isempty(tractset.analysispath)
    pth = fileparts(tractset.leadgroup);
    tractset.analysispath=[pth,filesep,'fiberfiltering',filesep,tractset.ID,'.fibfilt'];
end

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

% find which space and stim folder is used
if ~isfield(tractset.M,'pseudoM')
    stim_space = [ea_nt(tractset.native),'gs_',tractset.M.guid]; 
end

ea_dispercent(1/size(tractset.M.patient.list,1),'Computing E-field metrics on fibers')
for pt_i = 1:size(tractset.M.patient.list,1)

    %[~,subj_tag,~] = fileparts(tractset.M.patient.list{pt_i});
    %subSimPrefix = [subj_tag, '_sim-'];

    fprintf('\nProcessing: %s\n',tractset.M.patient.list{pt_i});

    for side = 1:2
    
        switch side
            case 1
                side_suffix = '_rh'; 
            case 2
                side_suffix = '_lh';
        end

        % OSS-DBS format (for Simbio, the function is not available)
        vatlist{pt_i,side};
        if strcmp(vatlist{pt_i,side},"skip")
            % no stimulation for this hemisphere
            continue
        elseif isfile(vatlist{pt_i,side})
            % for pseudoM, always re-define stim_space
            if isfield(tractset.M,'pseudoM')
                if tractset.M.pseudoM
                    [stim_folder_path,~,~] = fileparts(vatlist{pt_i,side});
                    [~,stim_folder,~] = fileparts(stim_folder_path);
                    stim_space = [ea_nt(tracetset.native),stim_folder];
                end
            end  

            % compute projection of the E-field onto the fibers
            ea_get_E_field_along_fibers(tractset.M.patient.list{pt_i}, stim_space, vatlist{pt_i,side}, merged_connectome, side_suffix, tractset.calcthreshold)
        else
            [~,pt_label,~] = fileparts(tractset.M.patient.list{pt_i});
            fprintf("Missing stimulation for %s, %s side \n",pt_label,side_suffix)
        end
    end

end
ea_dispercent(1/size(tractset.M.patient.list,1),'end')

