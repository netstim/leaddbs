function ea_getui(handles)

bids = getappdata(handles.leadfigure,'bids');
subjId = getappdata(handles.leadfigure,'subjId');

if ~isempty(bids)
    % Determine prefs path
    if isempty(handles.patientlist.Data) || isempty(handles.patientlist.Selection)
    	% prefsPath = fullfile(ea_getearoot, 'ea_ui.mat');
        prefsPath = '';
    else
    	prefsPath = bids.getPrefs(subjId{1}, 'uiprefs', 'mat');
    end

    if isfile(prefsPath)
        % Load UI prefs
        options = load(prefsPath);

        % Update UI
        ea_options2handles(options, handles);
    end
end
