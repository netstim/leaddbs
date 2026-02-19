function ea_init_connectome_tab(app_or_handles)
% Initialize Connectome tab in lead_dbs: parcellation list, FT methods,
% upsampling factor options, and default values (mirroring lead_connectome GUI).
% Call from lead_dbs_OpeningFcn with the app: ea_init_connectome_tab(app).
% Accepts either the lead_dbs app object (recommended) or handles struct.

% Use app component names when app is passed (has DropDown, FibertrackingMethodDropDown)
use_app = isobject(app_or_handles) && isprop(app_or_handles, 'DropDown');
if use_app
    app = app_or_handles;
    leadfig = app.leadfigure;
else
    handles = app_or_handles;
    if ~isfield(handles, 'parcellation') || ~isvalid(handles.parcellation)
        return
    end
    leadfig = handles.leadfigure;
end

earoot = ea_getearoot;
options.prefs = ea_prefs('');

%% 1. Parcellation scheme (same as lead_connectome OpeningFcn)
parcFiles = dir([ea_space([], 'labeling'), '*.nii']);
if ~isempty(parcFiles)
    parcellations = cellfun(@(x) {strrep(x, '.nii', '')}, {parcFiles.name});
    defaultParc = options.prefs.lc.defaultParcellation;
    parcIdx = find(ismember(parcellations, defaultParc), 1);
    if isempty(parcIdx), parcIdx = 1; end

    if use_app
        app.DropDown.Items = parcellations;
        app.DropDown.Value = parcellations{parcIdx};
    else
        if isprop(handles.parcellation, 'Items')
            handles.parcellation.Items = parcellations;
            handles.parcellation.Value = parcellations{parcIdx};
        else
            set(handles.parcellation, 'String', parcellations);
            set(handles.parcellation, 'Value', parcIdx);
        end
    end
end

%% 2. Fiber-tracking methods (same as lead_connectome: ea_ft_*.m + SPM version)
cnt = 1;
ndir = dir([earoot, 'connectomics', filesep, 'ea_ft_*.m']);
ftFunctions = cell(0);
ftMethods = cell(0);
for nd = length(ndir):-1:1
    [~, methodf] = fileparts(ndir(nd).name);
    try
        [thisndc, spmvers] = feval(methodf, 'prompt');
        if ismember(spm('ver'), spmvers)
            ftMethods{cnt} = thisndc;
            ftFunctions{cnt} = methodf;
            cnt = cnt + 1;
        end
    catch
        % skip methods that fail (e.g. wrong SPM)
    end
end
setappdata(leadfig, 'ftFunctions', ftFunctions);
if ~isempty(ftMethods)
    defaultFt = 'ea_ft_gqi_yeh';
    ftIdx = find(ismember(ftFunctions, defaultFt), 1);
    if isempty(ftIdx), ftIdx = 1; end
    if use_app
        app.FibertrackingMethodDropDown.Items = ftMethods;
        app.FibertrackingMethodDropDown.Value = ftMethods{ftIdx};
    else
        if isprop(handles.ftmethod, 'Items')
            handles.ftmethod.Items = ftMethods;
            try
                handles.ftmethod.Value = ftMethods{ftIdx};
            catch
                set(handles.ftmethod, 'Value', ftIdx);
            end
        else
            set(handles.ftmethod, 'String', ftMethods);
            set(handles.ftmethod, 'Value', ftIdx);
        end
    end
end

%% 3. Upsampling factor options (1..5 as in original)
upsamplingItems = {'1', '2', '3', '4', '5'};
if use_app
    app.UpsamplingFactorDropDown.Items = upsamplingItems;
    app.UpsamplingFactorDropDown.Value = '1';
else
    if isprop(handles.upsamplingfactor, 'Items')
        handles.upsamplingfactor.Items = upsamplingItems;
        handles.upsamplingfactor.Value = '1';
    else
        set(handles.upsamplingfactor, 'String', upsamplingItems);
        set(handles.upsamplingfactor, 'Value', 1);
    end
end

%% 4. Default values from prefs.machine.lc or ea_initlcopts
try
    lc = options.prefs.machine.lc;
catch
    lc = ea_initlcopts(leadfig);
end
if ~isfield(lc.struc.ft, 'upsample')
    lc.struc.ft.upsample.factor = 1;
    lc.struc.ft.upsample.how = 0;
end

% Apply defaults via handles (so lc2handles can set all checkboxes/edits)
if use_app
    [~, ~, handles] = ea_convertGUIDEArguments(app);
    ea_lc2handles(lc, handles);
    % Fiber Count visibility: show only for GQI (Yeh); set on app so it isn't lost (adapters may not support Visible)
    ftFunctions = getappdata(leadfig, 'ftFunctions');
    if ~isempty(ftFunctions)
        idx = find(ismember(app.FibertrackingMethodDropDown.Items, app.FibertrackingMethodDropDown.Value), 1);
        if isempty(idx), idx = 1; end
        showFiber = strcmp(ftFunctions{min(idx, length(ftFunctions))}, 'ea_ft_gqi_yeh');
        app.FiberCountEditField.Visible = showFiber;
        app.FiberCountLabel.Visible = showFiber;
    end
else
    ea_lc2handles(lc, handles);
end
