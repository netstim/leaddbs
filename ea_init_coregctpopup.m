function ea_init_coregctpopup(handles, defaultmethod, handlestring)

% Set coregctmethod popupmenu by default
% In checkreg window, we need to set coregmrmethod popupmenu
if ~exist('handlestring','var')
    handlestring = 'coregctmethod';
end

% Get functions and names
funcs = ea_regexpdir(ea_getearoot, 'ea_coregctmri_.*.m', 0);
funcs = regexp(funcs, '(ea_coregctmri_.*)(?=\.m)', 'match', 'once');
names = cellfun(@(x) eval([x, '(''prompt'');']), funcs, 'Uni', 0);

% Set names to popupmenu
set(handles.(handlestring), 'String', names);

% Set default method
set(handles.(handlestring), 'Value', find(ismember(names, defaultmethod), 1));
