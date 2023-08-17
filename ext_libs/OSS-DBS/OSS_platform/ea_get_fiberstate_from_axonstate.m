function ea_get_fiberstate_from_axonstate(ConnectomeName,ID_no_sub,stimdir,side,stimSetMode,native,convert_to_MNI)

    options.subj.subjID = ID_no_sub;
    options.native = native;
    %settings.connectome = ConnectomeName;
    settings.connectome = ['Multi-Tract: ', ConnectomeName];
    settings.stimSetMode = stimSetMode;


    %subDescPrefix = ['sub-', options.subj.subjId, '_desc-'];
    subSimPrefix = ['sub-', options.subj.subjID, '_sim-'];
    outputDir = stimdir;
    outputBasePath = [outputDir, filesep, subSimPrefix];
    % Extract connectome name
    connName = strrep(settings.connectome, 'Multi-Tract: ', '');
    settings.connectomePath = [outputDir, filesep, connName];



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

    axonState = ea_regexpdir([outputDir, filesep, 'Results_', sideCode], 'Axon_state.*\.mat', 0);
    if ~isempty(axonState)
        for f=1:length(axonState)

            % Determine tract name
            if startsWith(settings.connectome, 'Multi-Tract: ')
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
            ftr = load(axonState{f});
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

            % Save result for visualization

            % If stimSets, save to a corresponding folder
            if settings.stimSetMode
                resultProtocol = [outputDir, filesep, 'Result_StimProt_', sideStr, '_', stimProt_index];
                ea_mkdir(resultProtocol);
                if startsWith(settings.connectome, 'Multi-Tract: ')
                    fiberActivation = [resultProtocol, filesep, subSimPrefix, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '_tract-', tractName,'_prot-', stimProt_index, '.mat'];
                else
                    fiberActivation = [resultProtocol, filesep, subSimPrefix, 'fiberActivation_model-ossdbs_hemi-', sideLabel,'_prot-', stimProt_index, '.mat'];
                end
            else
                if startsWith(settings.connectome, 'Multi-Tract: ')
                    fiberActivation = [outputBasePath, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '_tract-', tractName, '.mat'];
                else
                    fiberActivation = [outputBasePath, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '.mat'];
                end
            end

            save(fiberActivation, '-struct', 'ftr');

            if convert_to_MNI==1 && options.native % Generate fiber activation file in MNI space
                fprintf('Restore connectome in MNI space: %s ...\n\n', settings.connectome);
                
                if startsWith(settings.connectome, 'Multi-Tract: ')
                    % load the particular pathway in MNI
                    [atlas_folder,~] = fileparts(tract);
                    originalFib = load([atlas_folder,filesep,tractName]);
                end
                conn = originalFib;

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
                    resultProtocol = [templateOutputDir, filesep, 'Result_StimProt_', sideStr, '_', stimProt_index];
                    ea_mkdir(resultProtocol);
                    if startsWith(settings.connectome, 'Multi-Tract: ')
                        fiberActivationMNI = [resultProtocol, filesep, subSimPrefix, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '_tract-', tractName,'_prot-', stimProt_index, '.mat'];
                    else
                        fiberActivationMNI = [resultProtocol, filesep, subSimPrefix, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '_prot-', stimProt_index, '.mat'];
                    end
                else
                    if startsWith(settings.connectome, 'Multi-Tract: ')
                        fiberActivationMNI = [templateOutputBasePath, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '_tract-', tractName, '.mat'];
                    else
                        fiberActivationMNI = [templateOutputBasePath, 'fiberActivation_model-ossdbs_hemi-', sideLabel, '.mat'];
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

        end
    end
end