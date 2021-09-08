function ea_storeui(handles)

if ~isfield(handles, 'patdir_choosebox')
    return;
end

bids = getappdata(handles.leadfigure,'bids');
subjId = getappdata(handles.leadfigure,'subjId');

% Determine prefs path
if strcmp(handles.patdir_choosebox.String, 'Choose Patient Directory') || length(subjId) > 1
	prefsPath = fullfile(ea_getearoot, 'ea_ui.mat');
else
	prefsPath = bids.getPrefs(subjId{1}, 'uiprefs', 'mat');
end

options = ea_handles2options(handles);

ea_mkdir(fileparts(prefsPath));
save(prefsPath, '-struct', 'options');
