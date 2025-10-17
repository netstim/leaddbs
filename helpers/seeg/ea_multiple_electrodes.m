function sel = ea_multiple_electrodes(varargin)

    models       = [];
    defaultCount = 2;     % safe fallback
    defaultModel = [];    % will fill from models later
    options = ea_handles2options(varargin{1});
    bids    = varargin{2};

    reco_file = fullfile( ...
        bids.datasetDir, 'derivatives','leaddbs', ...
        strcat('sub-', bids.subjId{1}), 'reconstruction', ...
        strcat('sub-', bids.subjId{1}, '_desc-reconstruction.mat'));

    % --- Try to load reco.props (elmodel, elname, labels) ---
    recoProps = try_load_reco_props(reco_file);  % [] if missing
    if ~isempty(recoProps)
        defaultCount = numel(recoProps);
        % If at least one elmodel is present, use that as the global default
        if isfield(recoProps, 'elmodel') && ~isempty(recoProps(1).elmodel)
            defaultModel = recoProps(1).elmodel;
        end
    end

    % --- Assemble model list (Lead-DBS if present; fallback hard-coded) ---
    models = try_get_models_from_leaddbs();
    models = normalize_model_list(models);
    if isempty(models)
        models = {'Medtronic 3389','Medtronic 3387','Abbott 6172 (Infinity)', ...
                  'Boston Sci Vercise','Boston Sci Cartesia (Directional)','Other'};
    end

    if isempty(defaultModel) || ~any(strcmp(models, defaultModel))
        defaultModel = models{1};
    end

    % Coerce count to a clean scalar int >=1
    defaultCount = max(1, round(coerce_double(defaultCount, 2)));

    sel = struct('count',[], 'models',{{}}, 'sameForAll',true, 'canceled',true);

    % ---------- Step 1: number of electrodes (wider dialog) ----------
    defStr = sprintf('%d', defaultCount);
    dlgOpts = struct('Resize','on','WindowStyle','normal','Interpreter','none');
    answ = inputdlg({'Number of electrodes:'}, 'Electrode Setup', [1 60], {defStr}, dlgOpts);
    if isempty(answ), return; end
    n = round(str2double(answ{1}));
    if ~isfinite(n) || n < 1, n = 1; end

    % ---------- Step 2: same model for all? ----------
    sameDefault = 'Yes';
    if ~isempty(recoProps)
        sameDefault = 'No';
    end
    ch = questdlg('Use the same electrode model for all electrodes?', ...
                  'Electrode Models','Yes','No','Cancel',sameDefault);
    if isempty(ch) || strcmpi(ch,'Cancel'), return; end
    same = strcmpi(ch,'Yes');

    % ---------- Step 3: pick model(s) ----------
    if same
        initIdx = find(strcmp(models, defaultModel),1);
        if isempty(initIdx), initIdx = 1; end

        idx = listdlg('PromptString','Select electrode model:', ...
                      'ListString',models, ...
                      'SelectionMode','single', ...
                      'InitialValue',initIdx, ...
                      'CancelString','Cancel', ...
                      'ListSize',[720 480]);   % wider/taller
        if isempty(idx), return; end

        sel.count      = n;
        sel.models     = repmat(models(idx), 1, n);
        sel.sameForAll = true;
        sel.canceled   = false;
    else
        % Per-electrode selection with informative multi-line prompt.
        per = cell(1,n);

        % Build per-electrode initial default indices if recoProps provided
        initIdxPer = ones(1,n);
        for i = 1:n
            if i <= numel(recoProps) && isfield(recoProps,'elmodel') && ~isempty(recoProps(i).elmodel)
                match = find(strcmp(models, aschar(recoProps(i).elmodel)), 1);
                if ~isempty(match)
                    initIdxPer(i) = match;
                else
                    initIdxPer(i) = 1;
                end
            else
                match = find(strcmp(models, defaultModel), 1);
                if isempty(match), match = 1; end
                initIdxPer(i) = match;
            end
        end

        for i = 1:n
            % Multi-line prompt: show Label and #contacts if available
            if i <= numel(recoProps)
                eln  = safe_elname(recoProps(i));
                nlab = safe_num_labels(recoProps(i));
                promptLines = { ...
                    sprintf('Model for electrode %d', i), ...
                    sprintf('Label: %s', eln), ...
                    sprintf('Number of contacts: %d', nlab) ...
                };
            else
                promptLines = {sprintf('Model for electrode %d', i)};
            end

            idx = listdlg('PromptString',promptLines, ...
                          'ListString',models, ...
                          'SelectionMode','single', ...
                          'InitialValue',initIdxPer(i), ...
                          'CancelString','Cancel', ...
                          'ListSize',[720 520]);   % wider/taller for clarity
            if isempty(idx), return; end
            per{i} = models{idx};
        end

        sel.count      = n;
        sel.models     = per;
        sel.sameForAll = false;
        sel.canceled   = false;
    end

    % ---------- Step 4: WRITE BACK into reco.props (if reco exists) ----------
    if exist(reco_file,'file') == 2
        try
            write_back_models_into_reco(reco_file, sel);
        catch ME
            warning('ea_multiple_electrodes:writeBackFailed', ...
                'Failed to write selected models back into reco.props: %s', ME.message);
        end
    end
end

% ================= helpers =================

function props = try_load_reco_props(matpath)
    props = [];
    try
        if exist(matpath,'file') == 2
            S = load(matpath);
            if isfield(S,'reco') && isstruct(S.reco) && isfield(S.reco,'props') && ~isempty(S.reco.props)
                raw = S.reco.props;
                % Normalize to struct array with fields: elmodel, elname, labels
                if istable(raw)
                    raw = table2struct(raw);
                end
                if isstruct(raw)
                    if iscolumn(raw), raw = raw(:)'; end
                    need = {'elmodel','elname','labels'};
                    for k = 1:numel(raw)
                        for f = 1:numel(need)
                            if ~isfield(raw, need{f}) || ~isfield(raw(k), need{f})
                                raw(k).(need{f}) = [];
                            end
                        end
                    end
                    props = raw;
                elseif iscell(raw)
                    tmp = struct('elmodel',[],'elname',[],'labels',[]);
                    P = repmat(tmp, 1, size(raw,1));
                    for k = 1:numel(P)
                        row = raw(k,:);
                        if numel(row) >= 1, P(k).elmodel = row{1}; end
                        if numel(row) >= 2, P(k).elname  = row{2}; end
                        if numel(row) >= 3, P(k).labels  = row{3}; end
                    end
                    props = P;
                end
            end
        end
    catch
        props = [];
    end
end

function write_back_models_into_reco(matpath, sel)
    % Load full MAT, update reco.props.elmodel (preserving props' original type), save.
    if exist(matpath,'file') ~= 2, return; end
    S = load(matpath);
    if ~isfield(S,'reco') || ~isstruct(S.reco)
        return; % nothing to do; user asked to only write if reco exists
    end
    if ~isfield(S.reco,'props') || isempty(S.reco.props)
        return; % no props to write back to
    end

    orig = S.reco.props;           % keep original type
    P    = to_struct_props(orig);  % struct array with elmodel/elname/labels

    % Resize P to sel.count
    n = sel.count;
    P = resize_props(P, n);

    % Assign selected models
    for i = 1:n
        P(i).elmodel = sel.models{i};
    end

    % Convert back to original type/shape
    S.reco.props = from_struct_props_like(P, orig);

    % Save the whole MAT back (preserve all variables in S)
    save(matpath, '-struct', 'S');
end

function P = to_struct_props(orig)
    % Return a struct array with fields elmodel, elname, labels
    if istable(orig)
        T = orig;
        % Expect columns named exactly 'elmodel','elname','labels'
        P = table2struct(T, 'ToScalar', false);
    elseif iscell(orig)
        tmp = struct('elmodel',[],'elname',[],'labels',[]);
        P = repmat(tmp, 1, size(orig,1));
        for k = 1:numel(P)
            row = orig(k,:);
            if ~isempty(row)
                if numel(row) >= 1, P(k).elmodel = row{1}; end
                if numel(row) >= 2, P(k).elname  = row{2}; end
                if numel(row) >= 3, P(k).labels  = row{3}; end
            end
        end
    elseif isstruct(orig)
        P = orig;
        if iscolumn(P), P = P(:)'; end
        need = {'elmodel','elname','labels'};
        for k = 1:numel(P)
            for f = 1:numel(need)
                if ~isfield(P, need{f}) || ~isfield(P(k), need{f})
                    P(k).(need{f}) = [];
                end
            end
        end
    else
        P = struct('elmodel',{},'elname',{},'labels',{});
    end
end

function P = resize_props(P, n)
    % Grow/shrink P to length n, preserving elname/labels where possible
    cur = numel(P);
    if cur == n, return; end
    if cur > n
        P = P(1:n);
    else
        template = struct('elmodel',[],'elname',[],'labels',[]);
        P(end+1:n) = template;
    end
end

function out = from_struct_props_like(P, like)
    % Convert struct array P back to the same representation as "like"
    if istable(like)
        % Make a table with columns elmodel, elname, labels
        % Ensure row vector -> column
        if isrow(P), P = P'; end
        elmodel = repmat({[]}, numel(P),1);
        elname  = repmat({[]}, numel(P),1);
        labels  = repmat({[]}, numel(P),1);
        for i = 1:numel(P)
            elmodel{i} = P(i).elmodel;
            elname{i}  = P(i).elname;
            labels{i}  = P(i).labels;
        end
        out = table(elmodel, elname, labels, 'VariableNames', {'elmodel','elname','labels'});
    elseif iscell(like)
        out = cell(numel(P), 3);
        for i = 1:numel(P)
            out{i,1} = P(i).elmodel;
            out{i,2} = P(i).elname;
            out{i,3} = P(i).labels;
        end
    elseif isstruct(like)
        % Keep as a row struct array
        out = P(:)';
    else
        % Fallback: write struct array
        out = P(:)';
    end
end

function name = safe_elname(p)
    name = 'unnamed';
    if isfield(p,'elname') && ~isempty(p.elname)
        if isstring(p.elname) || ischar(p.elname)
            name = char(p.elname);
        elseif iscell(p.elname) && ~isempty(p.elname) && (ischar(p.elname{1}) || isstring(p.elname{1}))
            name = char(p.elname{1});
        end
    end
end

function n = safe_num_labels(p)
    n = 0;
    if isfield(p,'labels') && ~isempty(p.labels)
        L = p.labels;
        if iscell(L)
            n = numel(L);
        elseif isstring(L)
            n = numel(L);
        elseif ischar(L)
            n = 1;
        elseif isnumeric(L)
            n = numel(L);
        end
    end
end

function tf = iscellstrlike(x)
    tf = iscell(x) && all(cellfun(@(c) ischar(c) || (isstring(c)&&isscalar(c)), x));
end

function tf = isnumericlike(x)
    tf = isnumeric(x) || (isstring(x)&&isscalar(x) && ~isnan(str2double(x))) ...
         || (ischar(x) && ~isempty(x) && ~isnan(str2double(x)));
end

function tf = istextlike(x)
    tf = ischar(x) || (isstring(x) && isscalar(x));
end

function y = aschar(x)
    if isstring(x), y = char(x); else, y = x; end
end

function y = coerce_double(x, fallback)
    try
        if isnumeric(x)
            y = double(x);
            if isempty(y) || ~isfinite(y), y = fallback; end
            if numel(y) > 1, y = y(1); end
        elseif isstring(x) || ischar(x)
            y = str2double(x);
            if ~isfinite(y), y = fallback; end
        else
            y = fallback;
        end
    catch
        y = fallback;
    end
end

function models = try_get_models_from_leaddbs()
    models = {};
    try
        if exist('ea_resolve_elspec','file')
            out = ea_resolve_elspec; % often 73x1 cell
            models = normalize_model_list(out);
        end
    catch
        % ignore
    end
end

function models = normalize_model_list(x)
    % Accepts cellstr (row/col), string array, char array
    % Returns row cellstr, unique, non-empty
    if isstring(x)
        x = cellstr(x);
    elseif ischar(x)
        x = cellstr(x);
    elseif ~iscell(x)
        x = {};
    end
    x = x(:); % column
    x = x(~cellfun(@(s) isempty(s) || all(isspace(char(s))), x)); % drop empty
    for k = 1:numel(x)
        if isstring(x{k}), x{k} = char(x{k}); end
    end
    if ~isempty(x)
        [~,ia] = unique(x, 'stable');
        x = x(sort(ia));
    end
    models = x(:)';  % force row
end