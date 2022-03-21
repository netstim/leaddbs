function ea_init_coregmrpopup(handles, defaultmethod)
% Initialize MR coregistration methods popupmenu.

% Set methods
names = {'ANTs (Avants 2008)'
    'BRAINSFit (Johnson 2007)'
    'FLIRT (Jenkinson 2001 & 2002)'
    'FLIRT BBR (Greve and Fischl 2009)'
    'Hybrid SPM & ANTs'
    'Hybrid SPM & BRAINSFIT'
    'Hybrid SPM & FLIRT'
    'SPM (Friston 2007)'};

prefs = ea_prefs;
if prefs.env.dev
    names = [names; 'ANTs Nonlinear Coregistration'];
end

% Set names to popupmenu
set(handles.coregmrmethod, 'String', names);

% Set default method
if ~ismember(defaultmethod, names)
    defaultmethod = 'SPM (Friston 2007)';
end

set(handles.coregmrmethod, 'Value', find(ismember(names, defaultmethod), 1));
