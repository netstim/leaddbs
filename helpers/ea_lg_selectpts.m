function ea_lg_selectpts(lgfile, selection, guid, outputfolder)
% Helper function to create a new Lead group analysis by selecting a subset
% of patients in the current group analysis and optionally set a new guid
% for the new group analysis

% Check M struct
if isstruct(lgfile)
    M = lgfile;
elseif isfile(lgfile)
    load(lgfile, 'M');
end

% Check selection
if ~exist('selection', 'var') || isempty(selection)
    selection = 1:length(M.patient.list);
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

% Select patients
M.patient.list = M.patient.list(selection);
M.patient.group = M.patient.group(selection);
if ~isempty(M.clinical.vars)
    for i=1:length(M.clinical.vars)
        M.clinical.vars{i} = M.clinical.vars{i}(selection, :);
    end
end
M.ui.listselect = M.ui.listselect(selection);
M.S = M.S(selection);
if ~isempty(M.isomatrix)
    for i=1:length(M.isomatrix)
        M.isomatrix{i} = M.isomatrix{i}(selection, :);
    end
end
if isfield(M, 'elstruct') && ~isempty(M.elstruct)
    M.elstruct = M.elstruct(selection);
end
if isfield(M, 'stats') && ~isempty(M.stats)
    M.stats = M.stats(selection);
end

% Fix group and group color
groupset = unique(M.patient.group);
M.groups.group = (1:length(groupset))';
M.groups.color = M.groups.color(groupset, :);
group = zeros(size(M.patient.group));
for i=1:length(groupset)
    group(M.patient.group==group(i)) = i;
end
M.patient.group = group;

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

% Copy stimulation folders
for i=1:length(M.patient.list)
    for nt=0:1
        stimFolder = [M.patient.list{i}, filesep, 'stimulations', filesep, ea_nt(nt)];
        if isfolder([stimFolder, 'gs_', oldGUID])
            if ~isfolder([stimFolder, 'gs_', guid])
                fprintf('Copying stimulation folder to:\n%s\n', [stimFolder, 'gs_', guid])
                copyfile([stimFolder, 'gs_', oldGUID], [stimFolder, 'gs_', guid]);
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
