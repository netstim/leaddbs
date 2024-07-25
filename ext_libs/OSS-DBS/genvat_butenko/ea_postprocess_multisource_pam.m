function ea_postprocess_multisource_pam(options,settings,side)

% Merge multisource PAMs using logical disjunction.
% IMPORTANT:fiberActivations here have the same pre-filtering with
% Phi_vector_max. That is why they contain the same fibers and can be
% merged directly
% By Butenko, konstantinmgtu@gmail.com

arguments
    options             % Lead-DBS options for electrode reconstruction and stimulation
    settings            % parameters for OSS-DBS simulation
    side                % 1 - rh, 2 - lh
end



% everything but tractname(end-5:end) should match when stacking!

% get all fiberActivations (make sure you don't have merged once before running it!)
fiberActivations = dir([settings.connectomeActivations,filesep,'sub-*']);

if side == 1
    % remove left side fiberActivations
    fiberActivations = fiberActivations(~ismember({fiberActivations.name}, {'-L.mat', '-L_'}));
else
    fiberActivations = fiberActivations(~ismember({fiberActivations.name}, {'-R.mat', '-R_'}));
end

already_merged_index = [];
for pam_idx = 1:size(fiberActivations,1)
    if ~ismember(pam_idx,already_merged_index)

        ftr = load([fiberActivations(pam_idx).folder,filesep,fiberActivations(pam_idx).name]);
        %[fibers_unique,ia,ic] = unique(ftr.fibers(:,4));

        fiberActivationsMerged = [fiberActivations(pam_idx).folder,filesep,fiberActivations(pam_idx).name(1:end-6),'.mat'];
        
        % search the rest to find the same pathway
        for pam_idx_2 = pam_idx:size(fiberActivations,1)
            % check if the same name except the source index
            if strcmp(fiberActivations(pam_idx).name(1:end-5),fiberActivations(pam_idx_2).name(1:end-5))

                % merge here
                ftr2 = load([fiberActivations(pam_idx_2).folder,filesep,fiberActivations(pam_idx_2).name]);
                % we only need to flip 0 to 1 if necessary. Other statutes
                % must stay the same per definition!

                ftr.fibers(ftr2.fibers(:,5) == 1,5) = 1;

                already_merged_index = [already_merged_index, pam_idx_2];

            end
        end
        % save the merged fiberActivation
        save(fiberActivationsMerged, '-struct', 'ftr');

        % for MNI, just load this particular fiberActivation and ovewrite
        % the column
        if options.native
            merged_status = ftr.fibers(:,5);
            ftr = load([fiberActivations(pam_idx).folder,filesep,fiberActivations(pam_idx).name]);
            ftr.fibers(:,5) = merged_status;
            fiberActivationsMergedMNI = [settings.connectomeActivationsMNI,filesep,fiberActivations(pam_idx).name(1:end-6),'.mat'];
            save(fiberActivationsMergedMNI, '-struct', 'ftr');
        end

    end

    % remove the source result to avoid wrong FF imports
    ea_delete([fiberActivations(pam_idx).folder,filesep,fiberActivations(pam_idx).name])

    if options.native
        if exist([settings.connectomeActivationsMNI,filesep,fiberActivations(pam_idx).name], 'file')
            ea_delete([settings.connectomeActivationsMNI,filesep,fiberActivations(pam_idx).name])
        end
    end

end

