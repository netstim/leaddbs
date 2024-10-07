function ea_postprocess_multisource_axonstates(axonStateFolder,options,settings)

% Merge multisource axon states using logical disjunction.
% IMPORTANT:this is a special case for multisource prob. PAM
% By Butenko, konstantinmgtu@gmail.com

arguments
    axonStateFolder     % folder where axons states across samples and sources are stored
    options             % Lead-DBS options for electrode reconstruction and stimulation
    settings            % parameters for OSS-DBS simulation
end


% everything but tractname(end-5:end) should match when stacking!

% get all fiberActivations (make sure you don't have merged once before running it!)
axonActivations = dir([axonStateFolder,filesep,'Axon_state_*']);

already_merged_index = [];
for pam_idx = 1:size(axonActivations,1)
    if ~ismember(pam_idx,already_merged_index)

        ftr = load([axonActivations(pam_idx).folder,filesep,axonActivations(pam_idx).name]);
        parts = strsplit(axonActivations(pam_idx).name,'_');
        source = parts{end}(1:end-4);
        sample = parts{end-1};
        % keep samples separate!
        axonActivationsMerged = [axonActivations(pam_idx).folder,filesep,axonActivations(pam_idx).name(1:end-8),'.mat'];

        % search the rest to find the same pathway from the same sample!
        for pam_idx_2 = pam_idx:size(axonActivations,1)
            % check if the same name and sample
            parts2 = strsplit(axonActivations(pam_idx_2).name,'_');
            source2 = parts2{end}(1:end-4);
            sample2 = parts2{end-1};

            if strcmp(axonActivations(pam_idx).name(1:end-8),axonActivations(pam_idx_2).name(1:end-8))

                % merge here
                ftr2 = load([axonActivations(pam_idx_2).folder,filesep,axonActivations(pam_idx_2).name]);
                % we only need to flip 0 to 1 if necessary. Other statutes
                % must stay the same per definition!
                ftr.fibers(ftr2.fibers(:,5) == 1,5) = 1;

                already_merged_index = [already_merged_index, pam_idx_2];

            end
        end
        % save the merged fiberActivation
        save(axonActivationsMerged, '-struct', 'ftr');

        % for MNI, just load this particular fiberActivation and ovewrite
        % the column
        % if options.native
        %     merged_status = ftr.fibers(:,5);
        %     ftr = load([fiberActivations(pam_idx).folder,filesep,fiberActivations(pam_idx).name]);
        %     ftr.fibers(:,5) = merged_status;
        %     fiberActivationsMergedMNI = [settings.connectomeActivationsMNI,filesep,fiberActivations(pam_idx).name(1:end-6),'.mat'];
        %     save(fiberActivationsMergedMNI, '-struct', 'ftr');
        % end

    end
    % 
    % remove the source result to avoid wrong FF imports
    ea_delete([axonActivations(pam_idx).folder,filesep,axonActivations(pam_idx).name])
    % 
    % if options.native
    %     if exist([settings.connectomeActivationsMNI,filesep,fiberActivations(pam_idx).name], 'file')
    %         ea_delete([settings.connectomeActivationsMNI,filesep,fiberActivations(pam_idx).name])
    %     end
    % end

end

