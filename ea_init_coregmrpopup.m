function ea_init_coregmrpopup(handles, defaultmethod)
% Initialize MR coregistration methods popupmenu.

% Get functions and names
funcs = ea_regexpdir(ea_getearoot, 'ea_coregpostopmr_.*.m', 0);
funcs = regexp(funcs, '(ea_coregpostopmr_.*)(?=\.m)', 'match', 'once');
names = cellfun(@(x) eval([x, '(''prompt'');']), funcs, 'Uni', 0);

% Set names to popupmenu
set(handles.coregmrmethod, 'String', names);

% Set default method
if ~ismember(defaultmethod, names)
    defaultmethod = 'SPM (Friston 2007)';
end

set(handles.coregmrmethod, 'Value', find(ismember(names, defaultmethod), 1));
