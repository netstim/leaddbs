function ea_lg_renewguid(lgfile, guid, outputfolder)
% Helper function to create a new Lead group analysis by setting a new guid
% to the input group analysis

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
    outputfolder = M.root;
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
save(ea_getGroupAnalysisFile(outputfolder), 'M');
