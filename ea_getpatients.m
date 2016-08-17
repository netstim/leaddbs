function uipatdir=ea_getpatients(handles)


p='/'; % default use root
try
p=pwd; % if possible use pwd instead (could not work if deployed)
end
try % finally use last patient parent dir if set.
earoot=ea_getearoot;
load([earoot,'ea_recentpatients.mat']);
p=fileparts(fullrpts{1});
end

uipatdir=ea_uigetdir(p,'Please choose patient folder(s)...');

if isempty(uipatdir)
    return
end

if exist('handles','var')
ea_load_pts(handles,uipatdir);
end