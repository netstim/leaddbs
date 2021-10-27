function ea_switch2develop

LeadRoot = ea_getearoot;

if isfile([LeadRoot, 'common', filesep, 'ea_recentgroups.mat'])
    disp('Backup recent groups ...');
    movefile([LeadRoot, 'common', filesep, 'ea_recentgroups.mat'], [LeadRoot, 'common', filesep, 'ea_recentgroups.mat.bids'])
end

if isfile([LeadRoot, 'common', filesep, 'ea_recentgroups.mat.dev'])
    disp('Restore recent groups ...');
    movefile([LeadRoot, 'common', filesep, 'ea_recentgroups.mat.dev'], [LeadRoot, 'common', filesep, 'ea_recentgroups.mat'])
end

if isfile([LeadRoot, 'common', filesep, 'ea_recentpatients.mat'])
    disp('Backup recent patients ...');
    movefile([LeadRoot, 'common', filesep, 'ea_recentpatients.mat'], [LeadRoot, 'common', filesep, 'ea_recentpatients.mat.bids'])
end

if isfile([LeadRoot, 'common', filesep, 'ea_recentpatients.mat.dev'])
    disp('Restore recent patients ...');
    movefile([LeadRoot, 'common', filesep, 'ea_recentpatients.mat.dev'], [LeadRoot, 'common', filesep, 'ea_recentpatients.mat'])
end

if isfile([LeadRoot, 'ea_ui.mat'])
    disp('Backup ea_ui.mat ...');
    movefile([LeadRoot, 'ea_ui.mat'], [LeadRoot, 'ea_ui.mat.bids'])
end

if isfile([LeadRoot, 'ea_ui.mat.dev'])
    disp('Restore ea_ui.mat ...');
    movefile([LeadRoot, 'ea_ui.mat.dev'], [LeadRoot, 'ea_ui.mat'])
end

if isfile([ea_gethome, '.ea_prefs.m'])
    disp('Backup .ea_prefs.m ...');
    movefile([ea_gethome, '.ea_prefs.m'], [ea_gethome, '.ea_prefs.m.bids'])
end

if isfile([ea_gethome, '.ea_prefs.m.dev'])
    disp('Restore .ea_prefs.m ...');
    movefile([ea_gethome, '.ea_prefs.m.dev'], [ea_gethome, '.ea_prefs.m'])
end

if isfile([ea_gethome, '.ea_prefs.mat'])
    disp('Backup .ea_prefs.mat ...');
    movefile([ea_gethome, '.ea_prefs.mat'], [ea_gethome, '.ea_prefs.mat.bids'])
end

if isfile([ea_gethome, '.ea_prefs.mat.dev'])
    disp('Restore .ea_prefs.mat ...');
    movefile([ea_gethome, '.ea_prefs.mat.dev'], [ea_gethome, '.ea_prefs.mat'])
end

if isfile([ea_gethome, '.ea_prefs.json'])
    disp('Backup .ea_prefs.json ...');
    movefile([ea_gethome, '.ea_prefs.json'], [ea_gethome, '.ea_prefs.json.bids'])
end

if isfile([ea_gethome, '.ea_prefs.json.dev'])
    disp('Restore .ea_prefs.json ...');
    movefile([ea_gethome, '.ea_prefs.json.dev'], [ea_gethome, '.ea_prefs.json'])
end
