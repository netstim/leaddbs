function ea_init_normpopup(handles, defaultmethod, handlestring)
% Initialize normalization methods popupmenu.

% Set normmethod popupmenu by default
% In checkreg window, we need to set coregmrmethod popupmenu
if ~exist('handlestring','var')
    handlestring = 'normmethod';
end

% Get functions and names
funcs = ea_regexpdir(ea_getearoot, 'ea_normalize_.*\.m$', 0);
funcs = regexp(funcs, '(ea_normalize_.*)(?=\.m)', 'match', 'once');
names = cellfun(@(x) eval([x, '(''prompt'');']), funcs, 'Uni', 0);

% Set names to popupmenu
set(handles.(handlestring), 'String', names);

% Set default method
if ~ismember(defaultmethod, names)
    defaultmethod = 'ANTs (Avants 2008)';
end

set(handles.(handlestring), 'Value', find(ismember(names, defaultmethod), 1));

% Check status of normalization setting button
ea_normsettings(handles, handlestring);
