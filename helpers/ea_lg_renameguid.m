function ea_lg_renameguid(lgfile, guid, outputfolder)
% Helper function to rename the input Lead group analysis by changing the
% guid. Stimulation folders will also be renamed.

% Check M struct
if isstruct(lgfile)
    M = lgfile;
elseif isfile(lgfile)
    load(lgfile, 'M');
end

oldGUID = M.guid;

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

% Adapt M.stats
if isfield(M, 'stats')
    for i=1:length(M.stats)
        for j=1:length(M.stats(i).ea_stats.stimulation.vat)
            M.stats(i).ea_stats.stimulation.vat(1).label = ['gs_', guid];
        end
        M.stats(1).ea_stats.stimulation.label = ['gs_', guid];
    end
end

% Move folders
for i=1:length(M.patient.list)
    for nt=0:1
        stimFolder = [M.patient.list{i}, filesep, 'stimulations', filesep, ea_nt(nt)];
        if isfolder([stimFolder, 'gs_', oldGUID])
            if ~isfolder([stimFolder, 'gs_', guid])
                fprintf('Moving stimulation folder to:\n%s\n', [stimFolder, 'gs_', guid])
                movefile([stimFolder, 'gs_', oldGUID], [stimFolder, 'gs_', guid]);
            else
                warning('off', 'backtrace');
                warning('Destination folder already exists! Skipping..\n%s', [stimFolder, 'gs_', guid]);
                warning('on', 'backtrace');
            end
        end
    end
end

% Save modified group analysis file
save(ea_getGroupAnalysisFile(outputfolder), 'M');
