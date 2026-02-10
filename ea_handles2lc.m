function lc=ea_handles2lc(handles)

% General settings (support both GUIDE and App Designer dropdowns)
parcellations = ea_lc_getdropdownitems(handles.parcellation);
parcVal = ea_lc_getdropdownvalue(handles.parcellation, parcellations);
lc.general.parcellation = parcellations{parcVal};

% Graph options:
lc.graph.struc_func_sim=get(handles.struc_func_sim,'Value');
lc.graph.nodal_efficiency=get(handles.nodal_efficiency,'Value');
lc.graph.eigenvector_centrality=get(handles.eigenvector_centrality,'Value');
lc.graph.degree_centrality=get(handles.degree_centrality,'Value');
lc.graph.fthresh=ea_lc_edit2num(handles.fthresh);
lc.graph.sthresh=ea_lc_edit2num(handles.sthresh);

% functional options:
lc.func.compute_CM=get(handles.compute_CM_func,'Value');
lc.func.compute_GM=get(handles.compute_GM_func,'Value');
lc.func.prefs.TR=ea_lc_edit2num(handles.TR);

% structural options:
lc.struc.compute_CM = get(handles.compute_CM_struc,'Value');
lc.struc.compute_GM = get(handles.compute_GM_struc,'Value');

% Fibertracking method list was stored on the figure
ftFunctions = getappdata(handles.leadfigure, 'ftFunctions');
if isempty(ftFunctions)
    ftFunctions = {}; % fallback
end
ftVal = ea_lc_getdropdownvalue(handles.ftmethod, ftFunctions);
if ~isempty(ftFunctions)
    lc.struc.ft.method = ftFunctions{ftVal};
else
    lc.struc.ft.method = '';
end

% Upsampling factor (dropdown may be missing in some layouts, default=1)
if isfield(handles, 'upsamplingfactor') && isvalid(handles.upsamplingfactor)
    lc.struc.ft.upsample.factor = ea_lc_getdropdownvalue(handles.upsamplingfactor, []);
else
    lc.struc.ft.upsample.factor = 1;
end

% Internal upsampling flag (checkbox optional, default=0)
if isfield(handles, 'use_internal_upsampling') && isvalid(handles.use_internal_upsampling)
    lc.struc.ft.upsample.how = get(handles.use_internal_upsampling,'Value');
else
    lc.struc.ft.upsample.how = 0;
end

lc.struc.ft.do        = get(handles.perf_ft,'Value');
lc.struc.ft.normalize = get(handles.normalize_fibers,'Value');

% Fiber count is only relevant / visible for certain FT methods
if isfield(handles,'fiber_count') && isvalid(handles.fiber_count) && ...
        strcmp(get(handles.fiber_count,'Visible'),'on')
    lc.struc.ft.dsistudio.fiber_count = ea_lc_edit2num(handles.fiber_count);
end

function items = ea_lc_getdropdownitems(h)
if isempty(h) || ~isvalid(h), items = {}; return; end
try
    if isprop(h, 'Items')
        items = h.Items;
        if isstring(items), items = cellstr(items); end
    else
        s = get(h, 'String');
        items = cellstr(s);
    end
catch
    items = {};
end

function idx = ea_lc_getdropdownvalue(h, items)
if isempty(h) || ~isvalid(h), idx = 1; return; end
try
    v = get(h, 'Value');
    if isnumeric(v) && v >= 1
        idx = double(v);
        return
    end
    if isempty(items)
        if isprop(h, 'Items'), items = h.Items; else items = get(h, 'String'); end
        if ischar(items), items = cellstr(items); end
    end
    if ischar(v) || isstring(v)
        idx = find(ismember(items, char(v)), 1);
        if isempty(idx), idx = 1; end
    else
        idx = 1;
    end
catch
    idx = 1;
end

function n = ea_lc_edit2num(h)
if isempty(h) || ~isvalid(h), n = nan; return; end
try
    if isprop(h, 'Value')
        n = h.Value;
    else
        n = str2double(get(h, 'String'));
    end
catch
    n = nan;
end
