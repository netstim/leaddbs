function ea_lg_adaptptsname(lgfile, ptsNames)
% Helper function to adapt patients name (folder) in Lead_groupanalysis.mat
% after patients' folder have been renamed

% Check M struct
if isstruct(lgfile)
    M = lgfile;
elseif isfile(lgfile)
    load(lgfile, 'M');
end

% Check new patient names
if length(M.patient.list) ~= ptsNames
    error('Length of the input name list doesn''t match the patient list in the group analysis file!');
end

% Adapat patient name: change from root/patient1A to root/patient1B
for i=1:length(M.patient.list)
    ptsRootFolder = fileparts(M.patient.list{i});
    M.patient.list{i} = fullfile(ptsRootFolder, ptsNames{i});
end

% Adapt M.elstruct
if isfield(M, 'elstruct')
    for i=1:length(M.elstruct)
        M.elstruct(i).name = ptsNames{i};
    end
end

% Adapt M.stats
if isfield(M, 'stats')
    for i=1:length(M.stats)
        for j=1:length(M.stats(i).ea_stats.patname)
            M.stats(i).ea_stats.patname{i} = ptsNames{i};
        end
        M.stats(i).ea_stats.electrodes.name = ptsNames{i};
    end
end

% Save modified group analysis file
save(ea_getGroupAnalysisFile(M.root), 'M');
