function ea_lg_adaptptsfolder(lgfile, ptsrootfolder)
% Helper function to adapt patients' root folder in Lead_groupanalysis.mat
% after moved to a new location

% Check M struct
if isstruct(lgfile)
    M = lgfile;
elseif isfile(lgfile)
    load(lgfile, 'M');
end

% Adapat patient folder: change from rootA/patient1 to rootB/patient1
if ~M.ui.detached % Avoid to do it for group analysis in detached mode
    for i=1:length(M.patient.list)
        [~, ptsName] = fileparts(M.patient.list{i});
        M.patient.list{i} = fullfile(ptsrootfolder, ptsName);
    end
end

% Save modified group analysis file
save(ea_getGroupAnalysisFile(M.ui.groupdir), 'M');
