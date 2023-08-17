function ea_normsettings(handles, handlestring)
% Set callback and Enable status of normalization setting button

% Set normmethod popupmenu by default
% In checkreg window, we need to set coregmrmethod popupmenu
if ~exist('handlestring','var')
    handlestring = 'normmethod';
end

% Get functions and names
funcs = ea_regexpdir(ea_getearoot, 'ea_normalize_.*\.m$', 0);
funcs = regexp(funcs, '(ea_normalize_.*)(?=\.m)', 'match', 'once');
names = cellfun(@(x) eval([x, '(''prompt'');']), funcs, 'Uni', 0);

% Get current normalization function
normMethod = handles.(handlestring).String{handles.(handlestring).Value};
normFunc = funcs{find(ismember(names, normMethod), 1)};

% Check if the method has setting available
[~, ~, hassettings] = feval(normFunc, 'prompt');

% Set setting button status
if hassettings
    set(handles.normsettings, 'Enable', 'on');
    % Set normsettings function name
    normsettingsfunc = strrep(normFunc, 'normalize', 'normsettings');
    setappdata(handles.normsettings, 'normsettingsfunc', normsettingsfunc);
else
    set(handles.normsettings, 'Enable', 'off');
end
