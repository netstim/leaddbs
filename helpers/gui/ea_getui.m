function ea_getui(handles)

bids = getappdata(handles.leadfigure,'bids');
subjId = getappdata(handles.leadfigure,'subjId');

% Determine prefs path
if strcmp(handles.patdir_choosebox.String, 'Choose Patient Directory') || length(subjId) > 1
	prefsPath = fullfile(ea_getearoot, 'ea_ui.mat');
else
	prefsPath = bids.getPrefs(subjId{1}, 'uiprefs', 'mat');
end

if isfile(prefsPath)
    % Load UI prefs
    options = load(prefsPath);

    % Update UI
    ea_options2handles(options, handles);
end
