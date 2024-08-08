function [settings,fibersFound] = ea_prepare_fibers(options, S, settings, outputPaths, source_i)
% Preprocess fibers: warp to native (if neccesary), remove those far away
% from the stimulating contacts and too short to allocate axons.
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    options     % Lead-DBS options for electrode reconstruction and stimulation
    S           % Lead-DBS stimulation settings
    settings    % parameters for OSS-DBS simulation
    outputPaths % various paths to conform with lead-dbs BIDS structure
    source_i            {mustBeNumeric} = 5; % source index. Not used if 5
end


% check if stimSets are used
if settings.stimSetMode
    if settings.optimizer || settings.trainANN
        stimProtocol{1,1} = string(ea_regexpdir([outputPaths.outputDir,filesep,'NB_rh'], '^Current_protocols_\d\.csv$', 0));
        stimProtocol{2,1} = string(ea_regexpdir([outputPaths.outputDir,filesep,'NB_lh'], '^Current_protocols_\d\.csv$', 0));
    else
        stimProtocol = ea_regexpdir(outputPaths.outputDir, '^Current_protocols_\d\.csv$', 0);
    end
else
    if ~settings.multisource
        stimProtocol = settings.Phi_vector;  % stim vector for source 1 only
    else
        stimProtocol = settings.Phi_vector_max;  % max stim vector accross sources
    end
end
%coords_mm = ea_load_reconstruction(options);
% load mni stimulation coordinates
load(options.subj.recon.recon, 'reco');
coords_mm_MNI = reco.('mni').coords_mm;

% path to a json with axon model description 
settings.pathwayParameterFile = 'Allocated_axons_parameters.json';

if isfield(options.subj, 'preopAnat')
    preopAnchor = options.subj.preopAnat.(options.subj.AnchorModality).coreg;
elseif options.native == 1 && ~isfield(options.subj, 'preopAnat')
    ea_warndlg('Native space info is missing, use template space instead')
    return
end

if ~startsWith(settings.connectome, 'Multi-Tract: ') % Normal connectome
    fprintf('Loading connectome: %s ...\n', settings.connectome);
    
    conn = load([ea_getconnectomebase, 'dMRI', filesep, settings.connectome, filesep, 'data.mat']);

    % Filter fibers based on the spherical ROI
    fiberFiltered = ea_filterfiber_stim(conn, coords_mm_MNI, stimProtocol, 'kuncel', 2);

    % Filter fibers based on the minimal length
    fiberFiltered = ea_filterfiber_len(fiberFiltered, settings.axonLength);

    % Move original fiber id to the 5th column, the 4th column will be 1:N
    fibersFound = zeros(size(fiberFiltered));
    for i=1:length(fiberFiltered)
        if ~isempty(fiberFiltered{i}.fibers)

            % Convert connectome fibers from MNI space to anchor space
            if options.native
                fprintf('Convert connectome into native space...\n\n');
                fibersMNIVox = ea_mm2vox(fiberFiltered{i}.fibers(:,1:3), [ea_space, options.primarytemplate, '.nii'])';
                fiberFiltered{i}.fibers(:,1:3)  = ea_map_coords(fibersMNIVox, ...
                    [ea_space, options.primarytemplate, '.nii'], ...
                    [options.subj.subjDir, filesep, 'forwardTransform'], ...
                    preopAnchor)';
            end

            fibers = zeros(size(fiberFiltered{i}.fibers,1),5);
            fibers(:,[1,2,3,5]) = fiberFiltered{i}.fibers;
            fibers(:,4) = repelem(1:length(fiberFiltered{i}.idx), fiberFiltered{i}.idx)';
            fiberFiltered{i}.fibers = fibers;
            fibersFound(i) = 1;
        end
    end

    settings.connectomePath = [outputPaths.outputDir, filesep, settings.connectome];
    ea_mkdir(settings.connectomePath);

    % also create a folder for PAM results
    settings.connectomeActivations = [settings.connectomePath,filesep,'PAM'];

    % clean up for the first source only
    if source_i == 1
        if exist(settings.connectomeActivations,'dir')
            ea_delete(settings.connectomeActivations);
        end
        ea_mkdir(settings.connectomeActivations);
    end

    if options.native
        settings.connectomePathMNI = [outputPaths.templateOutputDir, filesep, settings.connectome];
        if exist(settings.connectomePathMNI,'dir')
            ea_delete(settings.connectomePathMNI);
        end
        ea_mkdir(settings.connectomePathMNI);
        settings.connectomeActivationsMNI = [settings.connectomePathMNI,filesep,'PAM'];
        ea_mkdir(settings.connectomeActivationsMNI);
    end

    for i=1:length(fiberFiltered)
        % store the original number of fibers
        % to compute percent activation
        fiberFiltered{i}.origNum = size(conn.idx,1);
        buffer = fiberFiltered{i};
        save([settings.connectomePath, filesep, 'data', num2str(i), '.mat'], '-struct', 'buffer', '-v7.3');
    end
else % Multi-Tract connectome
    % Extract connectome name
    connName = strrep(settings.connectome, 'Multi-Tract: ', '');

    % Create output folder
    settings.connectomePath = [outputPaths.outputDir, filesep, connName];
    ea_mkdir(settings.connectomePath);

    % also create a folder for PAM results
    settings.connectomeActivations = [settings.connectomePath,filesep,'PAM'];
    % clean up for the first source only
    if source_i == 1
        if exist(settings.connectomeActivations,'dir')
            ea_delete(settings.connectomeActivations);
        end
        ea_mkdir(settings.connectomeActivations);

        if options.native
            settings.connectomePathMNI = [outputPaths.templateOutputDir, filesep, connName];
            if exist(settings.connectomePathMNI,'dir')
                ea_delete(settings.connectomePathMNI);
            end
            ea_mkdir(settings.connectomePathMNI);
            settings.connectomeActivationsMNI = [settings.connectomePathMNI,filesep,'PAM'];
            ea_mkdir(settings.connectomeActivationsMNI);
        end

    end

    % Get paths of tracts
    connFolder = [ea_getconnectomebase, 'dMRI_MultiTract', filesep, connName];
    tracts = ea_regexpdir(connFolder, '\.mat$', 0);

    settings.connectomeTractNames = cell(size(tracts));
    data1 = struct;
    data2 = struct;
    fibersFound = zeros(numel(tracts),2);
    for t=1:numel(tracts)
        tract = tracts{t};
        [~, tractName] = fileparts(tract);
        settings.connectomeTractNames{t} = tractName;
        fprintf('Loading connectome: %s, Tract: %s ...\n', connName, tractName);
        conn = load(tract);

        % Filter fibers based on the spherical ROI
        fiberFiltered = ea_filterfiber_stim(conn, coords_mm_MNI, stimProtocol, 'kuncel', 2);

        % Filter fibers based on the minimal length
        fiberFiltered = ea_filterfiber_len(fiberFiltered, settings.axonLength(t));

        % Move original fiber id to the 5th column, the 4th column will be 1:N
        for i=1:length(fiberFiltered)
            if ~isempty(fiberFiltered{i}.fibers)

                % Convert connectome fibers from MNI space to anchor space
                if options.native
                    fprintf('Convert connectome into native space...\n\n');
                    fibersMNIVox = ea_mm2vox(fiberFiltered{i}.fibers(:,1:3), [ea_space, options.primarytemplate, '.nii'])';
                    fiberFiltered{i}.fibers(:,1:3)  = ea_map_coords(fibersMNIVox, ...
                        [ea_space, options.primarytemplate, '.nii'], ...
                        [options.subj.subjDir, filesep, 'forwardTransform'], ...
                        preopAnchor)';
                end

                fibers = zeros(size(fiberFiltered{i}.fibers,1),5);
                fibers(:,[1,2,3,5]) = fiberFiltered{i}.fibers;
                fibers(:,4) = repelem(1:length(fiberFiltered{i}.idx), fiberFiltered{i}.idx)';
                fiberFiltered{i}.fibers = fibers;

                % store the original number of fibers
                % to compute percent activation
                fiberFiltered{i}.origNum = size(conn.idx,1);

                fibersFound(t,i) = 1;
            end

            if i == 1
                data1.(tractName) = fiberFiltered{1};
            else
                data2.(tractName) = fiberFiltered{2};
            end
            
        end
    end

    % Save filtered fibers
    save([settings.connectomePath, filesep, 'data1.mat'], '-struct', 'data1', '-v7.3');
    save([settings.connectomePath, filesep, 'data2.mat'], '-struct', 'data2', '-v7.3');
end
