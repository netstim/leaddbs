%% Working with new data and you do not want to use the GUI
% Setup the options struct.
options.root = fullfile(pwd);
options.patientname = 'P008';
options.uipatdirs = {fullfile(options.root, options.patientname)};
options.native = true;
options.loadnativereco = 0; % Load scrf reco when available
options.sides = [1 2];
options.prefs = ea_prefs;
options.elmodel = 'Medtronic 3389';
options = ea_resolve_elspec(options);

% Initialize the MERState, give it the options, and set it to a default state
temp1 = MERState();
temp1.setOptions(options);
temp1.clearData();
temp1.setDataToDefaults();

% Adjust for NexFrame rotation.
% We record three landmarks on the NexFrame in patient native space.
% These landmarks can then be used to adjust for the NexFrame rotation.
% These values are case-specific, do not use them in your own data!
nex_adj_left = struct('label', {'A', 'E', 'Entry'},...
    'coords', {[-50.3 72.8 117.0], [-49.6 50.4 128.1], [-34.8 36.6 72.5]});
nex_adj_right = struct('label', {'A', 'E', 'Entry'},...
    'coords', {[69.1 68.2 106.7], [68.8 46.3 118.2], [45.3 33.5 65.6]});
temp1.updateFrame('left', nex_adj_left);
temp1.updateFrame('right', nex_adj_right);

% Update into which MER track the DBS was implanted.
temp1.updateDBSDepth('left', 0.3);
temp1.updateDBSImplantTrack('left', 'medial');

% Add some markers.
temp1.addMarkerAtDepth('left', 'central', MERState.MarkerTypes.Top, '', 8.1);
temp1.addMarkerAtDepth('left', 'central', MERState.MarkerTypes.Bottom, '', 2.7);
temp1.addMarkerAtDepth('left', 'medial', MERState.MarkerTypes.Top, '', 6.5);
temp1.addMarkerAtDepth('left', 'medial', MERState.MarkerTypes.Bottom, '', 0.4);

% Save for later.
temp1.save('y');  % Pass 'y' to overwrite file if it exists. Omit to be prompted.

%% Working with data that were already saved (commandline or GUI)
temp2 = MERState();
temp2.Config.root = fullfile(pwd);
temp2.Config.patientname = 'P008';
temp2.load();
markers = temp2.exportMarkers()
