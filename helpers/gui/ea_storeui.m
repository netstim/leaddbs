function ea_storeui(handles)

if ~isfield(handles, 'patdir_choosebox') && ~isfield(handles, 'datasetselect')
    return;
end

bids = getappdata(handles.leadfigure,'bids');
subjId = getappdata(handles.leadfigure,'subjId');

% Determine prefs path
if length(subjId) > 1 || ...
        isfield(handles, 'datasetselect') && strcmp(handles.datasetselect.String, 'Choose Dataset Directory') || ...
        isfield(handles, 'patdir_choosebox') && strcmp(handles.patdir_choosebox.String, 'Choose Patient Directory')
	prefsPath = fullfile(ea_getearoot, 'ea_ui.mat');
else
	prefsPath = bids.getPrefs(subjId{1}, 'uiprefs', 'mat');
end

options = ea_handles2options(handles);

ea_mkdir(fileparts(prefsPath));
save(prefsPath, '-struct', 'options');
