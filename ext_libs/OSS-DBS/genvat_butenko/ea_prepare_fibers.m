function [settings,fibersFound] = ea_prepare_fibers(options, S, settings, outputPaths)

% check if classic S or stimSets are used
% note that PAM works only for one source
if settings.stimSetMode
    stimProtocol = ea_regexpdir(outputPaths.outputDir, '^Current_protocols_\d\.csv$', 0);
else
    stimProtocol = S;
end
coords_mm = ea_load_reconstruction(options);

preopAnchor = options.subj.preopAnat.(options.subj.AnchorModality).coreg;
if ~startsWith(settings.connectome, 'Multi-Tract: ') % Normal connectome
    fprintf('Loading connectome: %s ...\n', settings.connectome);
    conn = load([ea_getconnectomebase, 'dMRI', filesep, settings.connectome, filesep, 'data.mat']);
    if options.native
        %originalFib = conn;
        % Convert connectome fibers from MNI space to anchor space
        fprintf('Convert connectome into native space...\n\n');
        fibersMNIVox = ea_mm2vox(conn.fibers(:,1:3), [ea_space, options.primarytemplate, '.nii'])';
        conn.fibers(:,1:3)  = ea_map_coords(fibersMNIVox, ...
            [ea_space, options.primarytemplate, '.nii'], ...
            [options.subj.subjDir, filesep, 'forwardTransform'], ...
            preopAnchor)';
    end

    % Filter fibers based on the spherical ROI
    if options.native
	    fiberFiltered = ea_filterfiber_stim(conn, coords_mm, stimProtocol, 'kuncel', 2, preopAnchor);
    else
        fiberFiltered = ea_filterfiber_stim(conn, coords_mm, stimProtocol, 'kuncel', 2, [ea_space, options.primarytemplate, '.nii']);
    end

    % Filter fibers based on the minimal length
    fiberFiltered = ea_filterfiber_len(fiberFiltered, settings.axonLength);

    % Move original fiber id to the 5th column, the 4th column will be 1:N
    fibersFound = zeros(size(fiberFiltered));
    for i=1:length(fiberFiltered)
        if ~isempty(fiberFiltered{i}.fibers)
            fibers = zeros(size(fiberFiltered{i}.fibers,1),5);
            fibers(:,[1,2,3,5]) = fiberFiltered{i}.fibers;
            fibers(:,4) = repelem(1:length(fiberFiltered{i}.idx), fiberFiltered{i}.idx)';
            fiberFiltered{i}.fibers = fibers;
            fibersFound(i) = 1;
        end
    end

    settings.connectomePath = [outputPaths.outputDir, filesep, settings.connectome];
    ea_mkdir(settings.connectomePath);
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
        if options.native
            %originalFib = conn;
            % Convert connectome fibers from MNI space to anchor space
            fprintf('Convert connectome into native space...\n\n');
            fibersMNIVox = ea_mm2vox(conn.fibers(:,1:3), [ea_space, options.primarytemplate, '.nii'])';
            conn.fibers(:,1:3)  = ea_map_coords(fibersMNIVox, ...
                [ea_space, options.primarytemplate, '.nii'], ...
                [options.subj.subjDir, filesep, 'forwardTransform'], ...
                preopAnchor)';
        end

        % Filter fibers based on the spherical ROI
        if options.native
	        fiberFiltered = ea_filterfiber_stim(conn, coords_mm, stimProtocol, 'kuncel', 2, preopAnchor);
        else
            fiberFiltered = ea_filterfiber_stim(conn, coords_mm, stimProtocol, 'kuncel', 2, [ea_space, options.primarytemplate, '.nii']);
        end

        % Filter fibers based on the minimal length
        fiberFiltered = ea_filterfiber_len(fiberFiltered, settings.axonLength(t));

        % Move original fiber id to the 5th column, the 4th column will be 1:N
        for i=1:length(fiberFiltered)
            if ~isempty(fiberFiltered{i}.fibers)
                fibers = zeros(size(fiberFiltered{i}.fibers,1),5);
                fibers(:,[1,2,3,5]) = fiberFiltered{i}.fibers;
                fibers(:,4) = repelem(1:length(fiberFiltered{i}.idx), fiberFiltered{i}.idx)';
                fiberFiltered{i}.fibers = fibers;

                % store the original number of fibers
                % to compute percent activation
                fiberFiltered{i}.origNum = size(conn.idx,1);

                fibersFound(t,i) = 1;
            end
        end
        data1.(tractName) = fiberFiltered{1};
        data2.(tractName) = fiberFiltered{2};
    end

    % Save filtered fibers
    save([settings.connectomePath, filesep, 'data1.mat'], '-struct', 'data1', '-v7.3');
    save([settings.connectomePath, filesep, 'data2.mat'], '-struct', 'data2', '-v7.3');
end