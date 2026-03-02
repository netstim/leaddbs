function handles = ea_lc2handles(lc, handles)
% Apply Lead-Connectome options struct (lc) to GUI handles.
% Shared by lead_connectome and lead_dbs Connectome tab (ea_init_connectome_tab).

% General settings: parcellation
parcellations = ea_lc_getitems(handles.parcellation);
if ~isempty(parcellations)
    if ~isfield(lc.general, 'parcellation')
        options.prefs = ea_prefs('');
        defaultParc = options.prefs.lc.defaultParcellation;
        idx = find(ismember(parcellations, defaultParc), 1);
    else
        idx = find(ismember(parcellations, lc.general.parcellation), 1);
        if isempty(idx)
            options.prefs = ea_prefs('');
            idx = find(ismember(parcellations, options.prefs.lc.defaultParcellation), 1);
        end
    end
    if ~isempty(idx)
        ea_lc_setvalue(handles.parcellation, idx, parcellations);
    end
end

% Graph options
ea_lc_setcheck(handles, 'struc_func_sim', lc.graph.struc_func_sim);
ea_lc_setcheck(handles, 'nodal_efficiency', lc.graph.nodal_efficiency);
ea_lc_setcheck(handles, 'eigenvector_centrality', lc.graph.eigenvector_centrality);
ea_lc_setcheck(handles, 'degree_centrality', lc.graph.degree_centrality);
ea_lc_setedit(handles, 'fthresh', lc.graph.fthresh);
ea_lc_setedit(handles, 'sthresh', lc.graph.sthresh);

% Functional options
ea_lc_setcheck(handles, 'compute_CM_func', lc.func.compute_CM);
ea_lc_setcheck(handles, 'compute_GM_func', lc.func.compute_GM);
ea_lc_setedit(handles, 'TR', lc.func.prefs.TR);

% Structural options
ea_lc_setcheck(handles, 'compute_CM_struc', lc.struc.compute_CM);
ea_lc_setcheck(handles, 'compute_GM_struc', lc.struc.compute_GM);
ftFunctions = getappdata(handles.leadfigure, 'ftFunctions');
if ~isempty(ftFunctions)
    ftMethodIdx = find(ismember(ftFunctions, lc.struc.ft.method), 1);
    if isempty(ftMethodIdx)
        try
            defaultftMethod = ea_prefs('').machine.lc.struc.ft.method;
            ftMethodIdx = find(ismember(ftFunctions, defaultftMethod), 1);
        catch
            ftMethodIdx = 1;
        end
    end
    if isempty(ftMethodIdx), ftMethodIdx = 1; end
    ea_lc_setvalue(handles.ftmethod, ftMethodIdx, []);
end
if strcmp(lc.struc.ft.method, 'ea_ft_gqi_yeh')
    try set(handles.fiber_count, 'Visible', 'on'); end
    try set(handles.fiber_count_txt, 'Visible', 'on'); end
else
    try set(handles.fiber_count, 'Visible', 'off'); end
    try set(handles.fiber_count_txt, 'Visible', 'off'); end
end
ea_lc_setedit(handles, 'fiber_count', lc.struc.ft.dsistudio.fiber_count);
ea_lc_setcheck(handles, 'normalize_fibers', lc.struc.ft.normalize);
ea_lc_setcheck(handles, 'perf_ft', lc.struc.ft.do);
if isfield(lc.struc.ft, 'upsample')
    try
        ea_lc_setvalue(handles.upsamplingfactor, lc.struc.ft.upsample.factor, []);
    catch
    end
    ea_lc_setcheck(handles, 'use_internal_upsampling', lc.struc.ft.upsample.how);
end

function items = ea_lc_getitems(h)
if isempty(h) || ~isvalid(h)
    items = {};
    return
end
if isprop(h, 'Items')
    items = h.Items;
    if isstring(items), items = cellstr(items); end
else
    s = get(h, 'String');
    if ischar(s), s = {s}; end
    items = s;
end

function ea_lc_setvalue(h, idx, items)
if isempty(h) || ~isvalid(h), return; end
if isempty(items)
    items = ea_lc_getitems(h);
end
try
    if isprop(h, 'Value')
        if isempty(items)
            h.Value = idx;
        else
            h.Value = items{min(idx, length(items))};
        end
    else
        set(h, 'Value', idx);
    end
catch
end

function ea_lc_setcheck(handles, tag, val)
if ~isfield(handles, tag), return; end
try
    h = handles.(tag);
    if isvalid(h)
        if isprop(h, 'Value'), h.Value = val; else set(h, 'Value', val); end
    end
catch
end

function ea_lc_setedit(handles, tag, val)
if ~isfield(handles, tag), return; end
try
    h = handles.(tag);
    if isvalid(h)
        str = num2str(val);
        if isprop(h, 'Value')
            n = str2double(str);
            if ~isnan(n), h.Value = n; else h.Value = str; end
        else
            set(h, 'String', str);
        end
    end
catch
end
