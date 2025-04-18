function ea_unifiedmapping_getEproj(explorer,vatlist)
% Comptute E-field metrics directly on fibers (optionally in native space).
% The results are stored as
% patient_folder/connectomes/connnectomeName/stim_ID/E_metrics.mat (values for all fibers in the connectome!)
% By Butenko and Roediger, konstantinmgtu@gmail.com

arguments
    explorer       % fiber filtering object
    vatlist        % cell object {pt_N,2} with paths to 4D niftis
end

if isempty(explorer.analysispath)
    pth = fileparts(explorer.leadgroup);
    explorer.analysispath=[pth,filesep,'fiberfiltering',filesep,explorer.ID,'.fibfilt'];
end

% use connectome with all pathways combined
if explorer.multi_pathways == 1
    % check if merged_pathways is in fibfiltering folder
    [filepath,~,~] = fileparts(explorer.analysispath);
    merged_connectome = [filepath,filesep,explorer.calcsettings.connectome,filesep,'merged_pathways.mat'];
    if ~isfile(merged_connectome)
        % else check if it is in the original lead-group folder
        [filepath,~,~] = fileparts(explorer.leadgroup);
        merged_connectome = [filepath,filesep,explorer.calcsettings.connectome,filesep,'merged_pathways.mat'];
        if ~isfile(merged_connectome)
            % or if it is in another lead-group folder (where fibfiltering file is)
            [filepath,~,~] = fileparts(explorer.analysispath);
            [filepath,~,~] = fileparts(filepath);
            merged_connectome = [filepath,filesep,explorer.calcsettings.connectome,filesep,'merged_pathways.mat'];
        end
    end
else
    merged_connectome = [ea_getconnectomebase('dMRI'), explorer.calcsettings.connectome, filesep, 'data.mat'];
end

% find which space and stim folder is used
if ~isfield(explorer.M,'pseudoM')
    stim_space = [ea_nt(explorer.calcsettings.calcspace),'gs_',explorer.M.guid]; 
end

ea_dispercent(1/size(explorer.M.patient.list,1),'Computing E-field metrics on fibers')
for pt_i = 1:size(explorer.M.patient.list,1)

    %[~,subj_tag,~] = fileparts(tractset.M.patient.list{pt_i});
    %subSimPrefix = [subj_tag, '_sim-'];

    fprintf('\nProcessing: %s\n',explorer.M.patient.list{pt_i});

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
            if isfield(explorer.M,'pseudoM')
                if explorer.M.pseudoM
                    [stim_folder_path,~,~] = fileparts(vatlist{pt_i,side});
                    [~,stim_folder,~] = fileparts(stim_folder_path);
                    stim_space = [ea_nt(tracetset.calcsettings.calcspace),stim_folder];
                end
            end  

            % compute projection of the E-field onto the fibers
            ea_get_E_field_along_fibers(explorer.M.patient.list{pt_i}, stim_space, vatlist{pt_i,side}, merged_connectome, side_suffix, explorer.calcsettings.calcthreshold)
        else
            [~,pt_label,~] = fileparts(explorer.M.patient.list{pt_i});
            fprintf("Missing stimulation for %s, %s side \n",pt_label,side_suffix)
        end
    end

end
ea_dispercent(1/size(explorer.M.patient.list,1),'end')

