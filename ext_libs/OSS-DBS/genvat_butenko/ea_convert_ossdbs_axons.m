function ea_convert_ossdbs_axons(options,settings,side,prob_PAM,resultfig,outputPaths)
% Prepare Lead-DBS BIDS format fiber activations.
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    options             % Lead-DBS options for electrode reconstruction and stimulation
    settings            % parameters for OSS-DBS simulation
    side                {mustBeNumeric} % hemisphere index (0 - rh, 1 - lh)
    prob_PAM            {mustBeNumericOrLogical} % 1 if PAM is computed over an uncertain parameter (e.g. fiber diameter)
    resultfig           % figure handle
    outputPaths         % various paths to conform with lead-dbs BIDS structure 
end

switch side
    case 0
        sideLabel = 'R';
        sideCode = 'rh';
        sideStr = 'right';
    case 1
        sideLabel = 'L';
        sideCode = 'lh';
        sideStr = 'left';
end

if prob_PAM 
    axonStateProb = ea_regexpdir([outputPaths.outputDir, filesep, 'Results_', sideCode], 'Axon_state.*\_prob.mat', 0);
    axonState = axonStateProb;
else
    axonState = ea_regexpdir([outputPaths.outputDir, filesep, 'Results_', sideCode], 'Axon_state.*\.mat', 0);
end   

if ~isempty(axonState)
    for f=1:length(axonState)

        if prob_PAM
            axonState{f} = strrep(axonState{f},'_prob','');
        end

        % Determine tract name
        if startsWith(settings.connectome, 'Multi-Tract: ')
            connName = strrep(settings.connectome, 'Multi-Tract: ', '');
            tractName = regexp(axonState{f}, '(?<=Axon_state_).+(?=\.mat$)', 'match', 'once');
        end

        % If stimSetMode, extract the index from tractName (but axonState is still checked on the indexed file)
        if settings.stimSetMode
            if startsWith(settings.connectome, 'Multi-Tract: ')
                stimProt_index = regexp(tractName, '(?<=_)\d+$', 'match', 'once');
                tractName = regexp(tractName, '.+(?=_\d+$)', 'match', 'once');
            else
                stimProt_index = regexp(axonState{f}, '(?<=Axon_state_)\d+(?=\.mat$)', 'match', 'once');
            end
        end

        % Get fiber id and state from OSS-DBS result
        if prob_PAM
            ftr = load(axonStateProb{f});
        else
            ftr = load(axonState{f});
        end
        [fibId, ind] = unique(ftr.fibers(:,4));

        fibState = ftr.fibers(ind,5);

        % Restore full length fiber (as in original filtered fiber)
        if startsWith(settings.connectome, 'Multi-Tract: ')
            ftr = load([settings.connectomePath, filesep, 'data', num2str(side+1), '.mat'], tractName);
            ftr = ftr.(tractName);
            ftr.connectome_name = connName;
        else
            ftr = load([settings.connectomePath, filesep, 'data', num2str(side+1), '.mat']);
            ftr.connectome_name = settings.connectome;
        end
        ftr.fibers = ftr.fibers(ismember(ftr.fibers(:,4), fibId), :);
        originalFibID = ftr.fibers(:,5);

        % Extract original conn fiber id and idx, needed in case
        % calculation is done in native space
        [connFibID, idx] = unique(ftr.fibers(:,5));

        % Set fiber state
        for fib=1:length(fibId)
            ftr.fibers(ftr.fibers(:,4)==fibId(fib),5) = fibState(fib);
        end

        % Extract state of original conn fiber, needed in case
        % calculation is  done in native space
        connFibState = ftr.fibers(idx, 5);

        % Reset original fiber id as in the connectome
        ftr.fibers(:,4) = originalFibID;

        if ~prob_PAM
            % special status for prob_PAM is handled in ea_get_probab_axon_state() 
            if strcmp(settings.butenko_intersectStatus,'activated')
                ftr.fibers(ftr.fibers(:,5) == -1,5) = 1;
                ftr.fibers(ftr.fibers(:,5) == -3,5) = 1;
            elseif strcmp(settings.butenko_intersectStatus,'activated_at_active_contacts')
                ftr.fibers = OSS_DBS_Damaged2Activated(settings,ftr.fibers,ftr.idx,side+1);
            end
        end 

        % Save result for visualization

        % If stimSets, save to a corresponding folder
        if settings.stimSetMode
            resultProtocol = [outputPaths.outputDir, filesep, 'Result_StimProt_', sideStr, '_', stimProt_index];
            ea_mkdir(resultProtocol);
            if startsWith(settings.connectome, 'Multi-Tract: ')
                fiberActivation = [resultProtocol, filesep, subSimPrefix, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '_tract-', tractName,'_prot-', stimProt_index, '.mat'];
            else
                fiberActivation = [resultProtocol, filesep, subSimPrefix, 'fiberActivation_model-ossdbs_hemi-', sideLabel,'_prot-', stimProt_index, '.mat'];
            end
        else
            if startsWith(settings.connectome, 'Multi-Tract: ')
                fiberActivation = [outputPaths.outputBasePath, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '_tract-', tractName, '.mat'];
            else
                fiberActivation = [outputPaths.outputBasePath, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '.mat'];
            end
        end

        save(fiberActivation, '-struct', 'ftr');

        if options.native % Generate fiber activation file in MNI space
            fprintf('Restore connectome in MNI space: %s ...\n\n', settings.connectome);

            % if true_VTA, we should skip this step

            if startsWith(settings.connectome, 'Multi-Tract: ')
                % load the particular pathway in MNI
                %connName = strrep(settings.connectome, 'Multi-Tract: ', '');
                connFolder = [ea_getconnectomebase, 'dMRI_MultiTract', filesep, connName];
                conn = load([connFolder, filesep,tractName]);
            else
                conn = load([ea_getconnectomebase, 'dMRI', filesep, settings.connectome, filesep, 'data.mat']);
            end

            fprintf('Convert fiber activation result into MNI space...\n\n');
            conn.fibers = conn.fibers(ismember(conn.fibers(:,4), connFibID), :);
            % Set fiber state
            conn.fibers = [conn.fibers, zeros(size(conn.fibers,1),1)];
            for fib=1:length(connFibID)
                conn.fibers(conn.fibers(:,4)==connFibID(fib),5) = connFibState(fib);
            end

            % Recreate fiber idx
            [~, ~, idx] = unique(conn.fibers(:,4), 'stable');
            conn.idx = accumarray(idx,1);

            % Reset original fiber id as in the connectome
            ftr.fibers(:,4) = originalFibID;

            % Save MNI space fiber activation result

            % If stimSets, save to a corresponding folder
            if settings.stimSetMode
                resultProtocol = [outputPaths.templateOutputDir, filesep, 'Result_StimProt_', sideStr, '_', stimProt_index];
                ea_mkdir(resultProtocol);
                if startsWith(settings.connectome, 'Multi-Tract: ')
                    fiberActivationMNI = [resultProtocol, filesep, subSimPrefix, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '_tract-', tractName,'_prot-', stimProt_index, '.mat'];
                else
                    fiberActivationMNI = [resultProtocol, filesep, subSimPrefix, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '_prot-', stimProt_index, '.mat'];
                end
            else
                if startsWith(settings.connectome, 'Multi-Tract: ')
                    fiberActivationMNI = [outputPaths.templateOutputBasePath, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '_tract-', tractName, '.mat'];
                else
                    fiberActivationMNI = [outputPaths.templateOutputBasePath, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '.mat'];
                end
            end

            if startsWith(settings.connectome, 'Multi-Tract: ')
                conn.connectome_name = connName;
            else
                conn.connectome_name = settings.connectome;
            end

            save(fiberActivationMNI, '-struct', 'conn');

            if ~options.orignative % Visualize MNI space result
                fiberActivation = fiberActivationMNI;
            end
        end

        % Visualize fiber activation, but not for stimSetMode
        if ~settings.stimSetMode
            if ~isempty(resultfig)
                set(0, 'CurrentFigure', resultfig);
                if prob_PAM
                    ea_plot_prob_fiber_state(fiberActivation);
                else
                    ea_fiberactivation_viz(fiberActivation, resultfig);
                end
            end
        end
    end
end