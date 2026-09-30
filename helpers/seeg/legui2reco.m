
function legui2reco(options)
load(fullfile(options.root, options.patientname, 'reconstruction', ...
    strcat(options.patientname, '_electrodes.mat')));

%check
if isempty(ShaftMembership)
    error(['ShaftMembership is empty. Run AutoName, press Apply Changes, ' ...
           'close Assign Electrodes, and then save again.']);
end

if numel(ShaftMembership) ~= size(ElecXYZRaw,1)
    error(['ShaftMembership has %d entries, but there are %d contacts. ' ...
           'Run AutoName and Apply Changes again.'], ...
           numel(ShaftMembership), size(ElecXYZRaw,1));
end
% Parse labels and select depth electrodes
shaftIds = unique(ShaftMembership);
shaftIds = shaftIds(isfinite(shaftIds) & shaftIds > 0);
shaftIds = shaftIds(~isnan(shaftIds));

el_names = ShaftNames(shaftIds);

contact_counts = arrayfun(@(s) ...
    sum(ShaftMembership == s), shaftIds);
skip_mask = contact_counts < 4;

for idx = find(skip_mask)
    fprintf('Skipping electrode %s: only %d contacts found (<4); treating as artifact.\n', ...
        el_names{idx}, contact_counts(idx));
end

% Keep all shaft-related arrays aligned
keep_mask = ~skip_mask;

el_names       = el_names(keep_mask);
shaftIds       = shaftIds(keep_mask);
contact_counts = contact_counts(keep_mask);

if isempty(el_names)
    ea_warning('No electrodes with at least four contacts were found. Reconstruction was not updated.')
    return;
end

% --- Determine reconstruction path and check for existing file -----------
reco_dir      = fullfile(options.root, options.patientname, 'reconstruction');
ea_mkdir(reco_dir);
reco_filename = strcat(options.patientname, '_desc-reconstruction.mat');
reco_filepath = fullfile(reco_dir, reco_filename);

if exist(reco_filepath, 'file') == 2
    % Interactive prompt if in desktop MATLAB; default to "overwrite" in headless
    if usejava('desktop')
        answer = questdlg( ...
            sprintf('A reconstruction already exists for %s.\nRe-run automatic electrode detection and overwrite?', options.patientname), ...
            'Existing reconstruction found', ...
            'Yes (overwrite)','No (skip)','Cancel','No (skip)');
        if strcmpi(answer, 'Cancel')
            fprintf('Operation cancelled by user.\n');
            return;
        elseif strcmpi(answer, 'No (skip)')
            fprintf('Skipping automatic electrode assignment.\n');
            return;
        else
            fprintf('Overwriting existing reconstruction with new auto-electrode assignments...\n');
        end
    else
        % Headless: proceed
        fprintf('Existing reconstruction found; proceeding to overwrite (headless mode).\n');
    end
end

% --- Cache model list and pre-filter to SEEG candidates -------------------
allModels  = ea_resolve_elspec;            % list of model names
seegModels = filter_seeg_models(allModels);% exclude DBS; keep DIXI/AdTech/SEEG/etc.

% --- Build reco struct ----------------------------------------------------
reco = struct();
reco.props = struct('elmodel',{},'elname',{},'labels',{});
reco.native.coords_mm = cell(1, numel(el_names));
%reco.scrf.coords_mm   = cell(1, numel(el_names));
reco.mni.coords_mm    = cell(1, numel(el_names));
reco.native.markers   = struct('head',{},'tail',{},'x',{},'y',{});
%reco.scrf.markers     = struct('head',{},'tail',{},'x',{},'y',{});
reco.mni.markers      = struct('head',{},'tail',{},'x',{},'y',{});
reco.native.trajectory = repmat(struct(), 1, numel(el_names));
%reco.scrf.trajectory   = repmat(struct(), 1, numel(el_names));
reco.mni.trajectory    = repmat(struct(), 1, numel(el_names));

% --- Main loop over electrodes -------------------------------------------
for ii = 1:length(el_names)
    fprintf('Processing electrode %s.\n', el_names{ii})

    % Indices of contacts for this electrode
    shaftNum = shaftIds(ii);
    el_idx = find(ShaftMembership == shaftNum);
    if numel(el_idx) < 4
        warning('Skipping electrode %s because only %d contacts were found.', ...
            el_names{ii}, numel(el_idx));
        continue;
    end
    native_coords = ElecXYZRaw(el_idx, :);
%     proj_coords   = ElecXYZProjRaw(el_idx, :);
%     mni_coords    = ElecXYZMNIRaw(el_idx, :);

    % Decide start contact: lateral vs vertical
    if abs(max(native_coords(:,1)) - min(native_coords(:,1))) > 10
        [~, start_idx] = max(abs(native_coords(:,1)));
    else
        [~, start_idx] = max(native_coords(:,3));
    end
    d = sqrt(sum((native_coords - native_coords(start_idx,:)).^2, 2));
    [~, sort_idx] = sort(d, 'descend');

    % Fill coordinates
    reco.native.coords_mm{ii} = native_coords(sort_idx, :);
%     reco.scrf.coords_mm{ii}   = proj_coords(sort_idx, :);
%     reco.mni.coords_mm{ii}    = native_coords(sort_idx, :);

    % Label ordering sanity check
%     if ~issorted(c_numbers(el_idx(sort_idx)),'ascend')
%         ea_warning('Labels are inconsistent with automatic contact ordering.')
%     end

    % ======= AUTO-SELECT SEEG MODEL (by #contacts + mean spacing) =======
    n_contacts  = size(native_coords, 1);
    icd         = sqrt(sum(diff(native_coords(sort_idx,:)).^2, 2)); % inter-contact distances
    avg_spacing = mean(icd);                                        % mm

    % Use the model selected in the AutoName GUI
    shaftModelIdx = find(strcmp(ShaftNames, el_names{ii}), 1);

    if ~isempty(shaftModelIdx) && ...
            shaftModelIdx <= numel(ShaftModels) && ...
            ~isempty(ShaftModels{shaftModelIdx})

        elmodel = ShaftModels{shaftModelIdx};

        fprintf('Using selected model: %s for electrode %s\n', ...
            elmodel, el_names{ii});

    else
        % Fallback only if a saved model cannot be found
        elmodel = choose_best_seeg_model( ...
            seegModels, n_contacts, avg_spacing);

        fprintf(['No saved model found for %s. ' ...
            'Auto-selected model: %s\n'], ...
            el_names{ii}, elmodel);
    end

    % Resolve elspec for geometry-dependent calcs
    options.elmodel = elmodel;
    options = ea_resolve_elspec(options);

    % --- Write props
    reco.props(ii).elmodel             = elmodel;
    reco.props(ii).elname              = el_names{ii};
    reco.props(ii).labels              = ElecMapRaw((el_idx(sort_idx)),1);
    reco.props(ii).manually_corrected  = 1;

    % --- Markers & trajectories (native)
    reco.native.markers(ii).head = reco.native.coords_mm{ii}(1, :);
    reco.native.markers(ii).tail = reco.native.coords_mm{ii}(4, :);
    [xunitv, yunitv] = ea_calcxy_legui(reco.native.markers(ii).head, reco.native.markers(ii).tail);
    reco.native.markers(ii).x = reco.native.markers(ii).head + xunitv*(options.elspec.lead_diameter/2);
    reco.native.markers(ii).y = reco.native.markers(ii).head + yunitv*(options.elspec.lead_diameter/2);
    reco.native.trajectory(ii) = struct();
%     [~, reco.native.trajectory(ii), ~] = ea_resolvecoords(reco.native.markers(ii), elmodel);
%     reco.native.trajectory(ii) = ea_resolvecoords(reco.native.markers(ii), elmodel){1};
    % --- scrf
%     reco.scrf.markers(ii).head = reco.scrf.coords_mm{ii}(1, :);
%     reco.scrf.markers(ii).tail = reco.scrf.coords_mm{ii}(4, :);
%     [xunitv, yunitv] = ea_calcxy_legui(reco.scrf.markers(ii).head, reco.scrf.markers(ii).tail);
%     reco.scrf.markers(ii).x = reco.scrf.markers(ii).head + xunitv*(options.elspec.lead_diameter/2);
%     reco.scrf.markers(ii).y = reco.scrf.markers(ii).head + yunitv*(options.elspec.lead_diameter/2);
%     reco.scrf.trajectory(ii) = struct();
%     [~, reco.scrf.trajectory(ii), ~] = ea_resolvecoords(reco.scrf.markers(ii), elmodel);

%     % --- mni
%     reco.mni.markers(ii).head = reco.mni.coords_mm{ii}(1, :);
%     reco.mni.markers(ii).tail = reco.mni.coords_mm{ii}(4, :);
%     [xunitv, yunitv] = ea_calcxy_legui(reco.mni.markers(ii).head, reco.mni.markers(ii).tail);
%     reco.mni.markers(ii).x = reco.mni.markers(ii).head + xunitv*(options.elspec.lead_diameter/2);
%     reco.mni.markers(ii).y = reco.mni.markers(ii).head + yunitv*(options.elspec.lead_diameter/2);
%     try
%         reco.mni.trajectory(ii) = struct();
%         [~, reco.mni.trajectory(ii), ~] = ea_resolvecoords(reco.mni.markers(ii), elmodel);
%     catch
%         disp('Error building electrode trajectory');
%     end
end
% --- Apply approved Lead-DBS brain-shift transform -------------------------

disp(options.subj.brainshift)
brainshiftApproved = false;

if isfield(options, 'subj') && ...
        isfield(options.subj, 'brainshift') && ...
        isfield(options.subj.brainshift, 'log') && ...
        isfield(options.subj.brainshift.log, 'method') && ...
        isfile(options.subj.brainshift.log.method)

    bsjson = loadjson(options.subj.brainshift.log.method);

    if isfield(bsjson, 'approval') && bsjson.approval == 1
        brainshiftApproved = true;
    end
end

if brainshiftApproved && ...
        isfield(options.subj.brainshift, 'transform') && ...
        isfield(options.subj.brainshift.transform, 'scrf') && ...
        isfile(options.subj.brainshift.transform.scrf)

    d = load(options.subj.brainshift.transform.scrf);

    reco.scrf = applyScrfToLeGUIReco_local( ...
        d.mat, reco.native, 1:numel(el_names));

    fprintf('Applied APPROVED brain-shift transform to reconstruction.\n');

elseif brainshiftApproved
    warning(['Brain shift is approved, but the SCRF transform file ' ...
             'was not found. Brain shift was not applied.']);
else
    fprintf('No approved brain shift. Saving reconstruction without SCRF coordinates.\n');
end

if isfield(reco, 'scrf')
    sourceReco = reco.scrf;
    fprintf('Using SCRF-corrected coordinates for MNI normalization.\n');
else
    sourceReco = reco.native;
    fprintf('Using native coordinates for MNI normalization.\n');
end

% Transform reconstruction coordinates from anchorNative -> MNI


[whichnormmethod, template] = ...
    ea_whichnormmethod(options.subj.subjDir);

coregFile = ...
    options.subj.coreg.anat.preop.(options.subj.AnchorModality);

invDir = [options.subj.subjDir, filesep, 'inverseTransform'];

% Header for converting native world-mm coordinates to voxel coordinates
Vnative = spm_vol(coregFile);

reco.mni.coords_mm = cell(size(sourceReco.coords_mm));

for ii = 1:numel(sourceReco.coords_mm)

    coordsNativeMM = sourceReco.coords_mm{ii};

    if isempty(coordsNativeMM)
        continue
    end

    % native world-mm -> native voxel coordinates
    nativeVoxel = Vnative.mat \ ...
        [coordsNativeMM'; ones(1, size(coordsNativeMM,1))];

    nativeVoxel = nativeVoxel(1:3,:);

    % Transform ALL contacts for this electrode to MNI in one call
    mniCoords = ea_map_coords( ...
        nativeVoxel, ...
        coregFile, ...
        invDir, ...
        template, ...
        whichnormmethod);

    % ea_map_coords returns 3 x N; reco expects N x 3
    reco.mni.coords_mm{ii} = mniCoords';

end

%MNI Markers

reco.mni.markers = struct('head',{},'tail',{},'x',{},'y',{});

for ii = 1:numel(reco.mni.coords_mm)

    coords = reco.mni.coords_mm{ii};

    if isempty(coords)
        continue
    end

    reco.mni.markers(ii).head = coords(1,:);
    reco.mni.markers(ii).tail = coords(4,:);

    elmodel = reco.props(ii).elmodel;

    tmpOptions = options;
    tmpOptions.elmodel = elmodel;
    tmpOptions = ea_resolve_elspec(tmpOptions);

    [xunitv, yunitv] = ea_calcxy_legui( ...
        reco.mni.markers(ii).head, ...
        reco.mni.markers(ii).tail);

    reco.mni.markers(ii).x = ...
        reco.mni.markers(ii).head + ...
        xunitv*(tmpOptions.elspec.lead_diameter/2);

    reco.mni.markers(ii).y = ...
        reco.mni.markers(ii).head + ...
        yunitv*(tmpOptions.elspec.lead_diameter/2);
end

% --- Save reconstruction ---------------------------------------------------
save(reco_filepath, 'reco');

end % function


% ================= Helpers =================

function out = filter_seeg_models(allModels)
    if isstring(allModels), allModels = cellstr(allModels); end
    if ischar(allModels),   allModels = cellstr(allModels); end
    allModels = allModels(:);

    % Include SEEG vendors/keywords; exclude obvious DBS models
    inc = contains(lower(allModels), {'dixi','adtech','seeg','depth','pmg','strip','sde','grid'});
    exc = contains(lower(allModels), {'medtronic','boston','abbott','vercise','cartesia','3387','3389','infinity','directional'});
    out = allModels(inc & ~exc);

    if isempty(out)
        % fall back to everything but DBS if filter yields nothing
        out = allModels(~exc);
    end
    out = out(:)'; % row cellstr
end

function elmodel = choose_best_seeg_model(models, n_contacts_obs, spacing_obs)
    % Score candidates by contact count and spacing similarity.
    bestScore = inf;
    bestModel = models{1};

    for k = 1:numel(models)
        mdl = models{k};

        % Query elspec for this model
        tmp.elmodel = mdl;
        try
            tmp = ea_resolve_elspec(tmp);
            spec = [];
            if isfield(tmp, 'elspec'), spec = tmp.elspec; end
        catch
            spec = [];
        end

        n_pred   = get_first_numeric(spec, {'n_contacts','ncontacts','contacts','num_contacts','numelectrodes'});
        pitch_mm = get_first_numeric(spec, {'spacing','intercontact','inter_contact_distance','ringpitch','contact_distance','contact_pitch'});

        % Contact-count penalty
        if ~isempty(n_pred) && isfinite(n_pred)
            p_count = abs(double(n_pred) - double(n_contacts_obs));
        else
            % infer from model name if possible
            nums = regexp(mdl,'\d+','match');
            if ~isempty(nums)
                n_infer = str2double(nums{end});
                if isfinite(n_infer)
                    p_count = abs(n_infer - double(n_contacts_obs));
                else
                    p_count = 2; % mild penalty
                end
            else
                p_count = 2; % mild penalty
            end
        end

        % Spacing penalty (mm)
        if ~isempty(pitch_mm) && isfinite(pitch_mm) && isfinite(spacing_obs)
            p_space = abs(double(pitch_mm) - double(spacing_obs));
        else
            p_space = 1.0; % neutral penalty if unknown
        end

        % Combined score
        score = 1.5*p_count + 1.0*p_space;

        % Prefer SEEG vendors on tie
        isPreferred = contains(lower(mdl), {'dixi','adtech','seeg'});
        if score < bestScore || (abs(score-bestScore) < 1e-6 && isPreferred)
            bestScore = score;
            bestModel = mdl;
        end
    end

    elmodel = bestModel;
end

function val = get_first_numeric(spec, names)
    val = [];
    if isempty(spec) || ~isstruct(spec), return; end
    for i = 1:numel(names)
        f = names{i};
        if isfield(spec, f)
            v = spec.(f);
            if isnumeric(v) && isscalar(v) && isfinite(v)
                val = double(v);
                return;
            end
        end
    end
end

function scrf = applyScrfToLeGUIReco_local(mat, native, shafts)

scrf = native;

for shaft = shafts

    % Coordinates
    coords = native.coords_mm{shaft};

    transformed = mat * ...
        [coords, ones(size(coords,1),1)]';

    scrf.coords_mm{shaft} = transformed(1:3,:)';

    % Markers
    markerNames = {'head','tail','x','y'};

    for k = 1:numel(markerNames)

        fieldName = markerNames{k};
        point = native.markers(shaft).(fieldName);

        if isempty(point)
            scrf.markers(shaft).(fieldName) = [];
            continue;
        end

        transformed = mat * [point, 1]';

        scrf.markers(shaft).(fieldName) = ...
            transformed(1:3)';
    end
end
end
