function startPath = ea_startpath
% Get start path for uigetfile/uigetdir

startPath = ea_gethome; % Use home folder by default

try
    startPath = pwd; % if possible use pwd instead (could not work if deployed)
end

try % Finally try to use the parent folder of the last recent patient folder
    load([ea_getearoot, 'common', filesep, 'ea_recentpatients.mat'], 'recentfolders');
    startPath = fileparts(recentfolders{1});
end
