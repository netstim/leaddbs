function p = ea_startpath
p='/'; % default use root
try
    p=pwd; % if possible use pwd instead (could not work if deployed)
end
try % finally use last patient parent dir if set.
    load([ea_prefsdir, filesep, 'ea_recentpatients.mat']);
    p=fileparts(fullrpts{1});
end