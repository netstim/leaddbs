function ea_lg_setguid(lgfile, guid, outputfolder)
% Helper function to set a new guid for Lead group analysis

% Check M struct
if isstruct(lgfile)
    M = lgfile;
elseif isfile(lgfile)
    load(lgfile, 'M');
end

% Check guid
if ~exist('guid', 'var') || isempty(guid)
    guid = datestr(datevec(now), 'yyyymmddHHMMSS');
end

% Check output folder
if ~exist('outputfolder', 'var') || isempty(outputfolder)
    outputfolder = M.ui.groupdir;
else
    if ~isfolder(outputfolder)
        mkdir(outputfolder);
    end
end

% Set guid
M.guid = guid;
if ~isempty(M.S)
    for i=1:length(M.S)
        M.S(i).label = ['gs_', guid];
    end
end

% Remove stats for the new group analysis file.
if isfield(M, 'stats')
    M = rmfield(M, 'stats');
end

% Save modified group analysis file
save(fullfile(outputfolder, 'LEAD_groupanalysis.mat'), 'M');
