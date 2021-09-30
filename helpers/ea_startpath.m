function p = ea_startpath
p='/'; % default use root
try
    p=pwd; % if possible use pwd instead (could not work if deployed)
end
try % finally use last patient parent dir if set.
    earoot=ea_getearoot;
    load([earoot,'common',filesep,'ea_recentpatients.mat']);
    p=fileparts(recentfolders{1});
end
