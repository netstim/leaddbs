
function legui2reco(options)
load(fullfile(options.root, options.patientname, 'reconstruction', ...
     strcat(options.patientname, '_electrodes.mat')));

% Parse labels and select depth electrodes
c_labels  = cellfun(@(x) regexprep(x, '\d+$', ''), ElecMapRaw(:,1), 'UniformOutput', false);
c_numbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match')), ElecMapRaw(:,1));
el_names  = unique(c_labels(DepthElecRaw)); % only depth electrodes (SEEG)
contact_counts = cellfun(@(name) sum(strcmp(c_labels, name)), el_names);
skip_mask = contact_counts < 4;

for idx = find(skip_mask)
    fprintf('Skipping electrode %s: only %d contacts found (<4); treating as artifact.\n', ...
        el_names{idx}, contact_counts(idx));
end

el_names = el_names(~skip_mask);

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
reco.scrf.coords_mm   = cell(1, numel(el_names));
reco.mni.coords_mm    = cell(1, numel(el_names));
reco.native.markers   = struct('head',{},'tail',{},'x',{},'y',{});
reco.scrf.markers     = struct('head',{},'tail',{},'x',{},'y',{});
reco.mni.markers      = struct('head',{},'tail',{},'x',{},'y',{});
reco.native.trajectory = repmat(struct(), 1, numel(el_names));
reco.scrf.trajectory   = repmat(struct(), 1, numel(el_names));
reco.mni.trajectory    = repmat(struct(), 1, numel(el_names));

% --- Main loop over electrodes -------------------------------------------
for ii = 1:length(el_names)
    fprintf('Processing electrode %s.\n', el_names{ii})

    % Indices of contacts for this electrode
    el_idx        = find(strcmp(c_labels, el_names{ii}));
    native_coords = ElecXYZRaw(el_idx, :);
    proj_coords   = ElecXYZProjRaw(el_idx, :);
    mni_coords    = ElecXYZMNIRaw(el_idx, :);

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
    reco.scrf.coords_mm{ii}   = proj_coords(sort_idx, :);
    reco.mni.coords_mm{ii}    = native_coords(sort_idx, :);

    % Label ordering sanity check
    if ~issorted(c_numbers(el_idx(sort_idx)),'ascend')
        ea_warning('Labels are inconsistent with automatic contact ordering.')
    end

    % ======= AUTO-SELECT SEEG MODEL (by #contacts + mean spacing) =======
    n_contacts  = size(native_coords, 1);
    icd         = sqrt(sum(diff(native_coords(sort_idx,:)).^2, 2)); % inter-contact distances
    avg_spacing = mean(icd);                                        % mm

    elmodel = choose_best_seeg_model(seegModels, n_contacts, avg_spacing);
    fprintf('Auto-selected model: %s (n=%d, mean spacing=%.2f mm)\n', elmodel, n_contacts, avg_spacing);

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
    [xunitv, yunitv] = ea_calcxy(reco.native.markers(ii).head, reco.native.markers(ii).tail);
    reco.native.markers(ii).x = reco.native.markers(ii).head + xunitv*(options.elspec.lead_diameter/2);
    reco.native.markers(ii).y = reco.native.markers(ii).head + yunitv*(options.elspec.lead_diameter/2);
    reco.native.trajectory(ii) = struct();
%     [~, reco.native.trajectory(ii), ~] = ea_resolvecoords(reco.native.markers(ii), elmodel);
%     reco.native.trajectory(ii) = ea_resolvecoords(reco.native.markers(ii), elmodel){1};
    % --- scrf
    reco.scrf.markers(ii).head = reco.scrf.coords_mm{ii}(1, :);
    reco.scrf.markers(ii).tail = reco.scrf.coords_mm{ii}(4, :);
    [xunitv, yunitv] = ea_calcxy(reco.scrf.markers(ii).head, reco.scrf.markers(ii).tail);
    reco.scrf.markers(ii).x = reco.scrf.markers(ii).head + xunitv*(options.elspec.lead_diameter/2);
    reco.scrf.markers(ii).y = reco.scrf.markers(ii).head + yunitv*(options.elspec.lead_diameter/2);
%     reco.scrf.trajectory(ii) = struct();
%     [~, reco.scrf.trajectory(ii), ~] = ea_resolvecoords(reco.scrf.markers(ii), elmodel);

    % --- mni
    reco.mni.markers(ii).head = reco.mni.coords_mm{ii}(1, :);
    reco.mni.markers(ii).tail = reco.mni.coords_mm{ii}(4, :);
    [xunitv, yunitv] = ea_calcxy(reco.mni.markers(ii).head, reco.mni.markers(ii).tail);
    reco.mni.markers(ii).x = reco.mni.markers(ii).head + xunitv*(options.elspec.lead_diameter/2);
    reco.mni.markers(ii).y = reco.mni.markers(ii).head + yunitv*(options.elspec.lead_diameter/2);
%     try
%         reco.mni.trajectory(ii) = struct();
%         [~, reco.mni.trajectory(ii), ~] = ea_resolvecoords(reco.mni.markers(ii), elmodel);
%     catch
%         disp('Error building electrode trajectory');
%     end
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
