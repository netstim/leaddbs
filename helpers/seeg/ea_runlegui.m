function run = ea_runlegui(options)
    addpath(fullfile(options.earoot, 'ext_libs', 'LeGUI', 'icons'));
	prefsPath = options.bids.getPrefs(options.bids.subjId{1}, 'uiprefs', 'mat');
    prefs = load(prefsPath);
    prefs.reconmethod = 'LeGUI (Davis 2021)';
    save(prefsPath, '-struct', 'prefs');
    app = LeGUI(fullfile(options.root, options.patientname), options);
end