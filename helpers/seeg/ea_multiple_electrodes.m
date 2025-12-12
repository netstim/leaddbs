function sel = ea_multiple_electrodes(varargin)
% EA_MULTIPLE_ELECTRODES
% Standalone popup(s) to pick:
%   1) Number of electrodes
%   2) Electrode model (one-for-all OR per-electrode)
%
% RETURNS struct SEL:
%   .count      double              % number of electrodes
%   .models     1 x count cellstr   % model per electrode
%   .sameForAll logical
%   .canceled   logical
%
% CALL STYLES (all optional; super forgiving):
%   sel = ea_multiple_electrodes();                          % auto get models from ea_resolve_elspec
%   sel = ea_multiple_electrodes(models);                    % provide model list (cellstr)
%   sel = ea_multiple_electrodes(models, defaultCount);      % also set default #electrodes
%   sel = ea_multiple_electrodes(models, defaultCount, defaultModel);
%
% It also tolerates calls like ea_multiple_electrodes(handles, bids) and will
% just ignore those and use defaults.

    % ---------- parse forgiving inputs ----------
    models       = [];
    defaultCount = 2;     % safe fallback
    defaultModel = [];    % will fill from models later

    % If the first arg looks like a model list, use it. Otherwise ignore.
    if nargin >= 1 && iscellstrlike(varargin{1})
        models = varargin{1};
    end

    % If second arg looks numeric-like, use it as defaultCount
    if nargin >= 2 && isnumericlike(varargin{2})
        defaultCount = coerce_double(varargin{2}, 2);
    end

    % Third arg as defaultModel if text-like
    if nargin >= 3 && istextlike(varargin{3})
        defaultModel = aschar(varargin{3});
    end

    % ---------- resolve/normalize model list ----------
    if isempty(models)
        models = try_get_models_from_leaddbs();
    end
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

    % ---------- Step 1: number of electrodes ----------
    defStr = sprintf('%d', defaultCount);     % <- safe even if defaultCount came in odd
    answ = inputdlg({'Number of electrodes:'}, 'Electrode Setup', 1, {defStr});
    if isempty(answ), return; end
    n = round(str2double(answ{1}));
    if ~isfinite(n) || n < 1, n = 1; end

    % ---------- Step 2: same model for all? ----------
    ch = questdlg('Use the same electrode model for all electrodes?', ...
                  'Electrode Models','Yes','No','Cancel','Yes');
    if isempty(ch) || strcmpi(ch,'Cancel'), return; end
    same = strcmpi(ch,'Yes');

    initIdx = find(strcmp(models, defaultModel),1);
    if isempty(initIdx), initIdx = 1; end

    % ---------- Step 3: pick model(s) ----------
    if same
        idx = listdlg('PromptString','Select electrode model:', ...
                      'ListString',models,'SelectionMode','single', ...
                      'InitialValue',initIdx,'CancelString','Cancel');
        if isempty(idx), return; end
        sel.count      = n;
        sel.models     = repmat(models(idx), 1, n);
        sel.sameForAll = true;
        sel.canceled   = false;
    else
        per = cell(1,n);
        for i = 1:n
            idx = listdlg('PromptString',sprintf('Model for electrode %d:', i), ...
                          'ListString',models,'SelectionMode','single', ...
                          'InitialValue',initIdx,'CancelString','Cancel');
            if isempty(idx), return; end
            per{i} = models{idx};
        end
        sel.count      = n;
        sel.models     = per;
        sel.sameForAll = false;
        sel.canceled   = false;
    end
end

% ================= helpers =================

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
    
    % Flatten into column
    x = x(:);
    
    % Remove empty / whitespace-only
    x = x(~cellfun(@(s) isempty(s) || all(isspace(char(s))), x));
    
    % Ensure everything is char (not string)
    for k = 1:numel(x)
        if isstring(x{k}), x{k} = char(x{k}); end
    end
    
    % Keep unique values in stable order
    if ~isempty(x)
        [~,ia] = unique(x, 'stable');
        x = x(sort(ia));
    end
    
    % Force row
    models = x(:)';  
end