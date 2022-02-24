function [apps, ids, accs] = ea_getsuiteapps(prefs)

if ~exist('prefs','var')
    prefs = ea_prefs;
end

apps = {'Lead DBS','Lead Connectome','Lead Anatomy','Lead Connectome Mapper','Lead Group'};
ids = {'lead_dbs','lead_connectome','lead_anatomy','lead_mapper','lead_group'};
accs = {'D','C','T','M','G'};

if prefs.env.dev
    apps = [apps, 'Lead Group Connectome','Lead OR','Lead Predict'];
    ids = [ids, 'lead_group_connectome','lead_or','lead_predict'];
    accs = [accs, 'K','O','I'];
end
