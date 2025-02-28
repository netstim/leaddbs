function ea_rebasegrouppts(handles)
% Adapt patient base folder ('dataset/derivatives/leaddbs' folder)

newBaseFolder = uigetdir('', 'Please select new base folder...');

M = getappdata(handles.leadfigure, 'M');
for i=1:length(M.patient.list)
    [~, ptsName] = fileparts(M.patient.list{i});
    M.patient.list{i} = fullfile(newBaseFolder, ptsName);
end
setappdata(handles.leadfigure, 'M', M);

ea_refresh_lg(handles);

analysisFile = ea_regexpdir(M.root, ['^dataset-[^\W_]+_analysis-', M.guid, '\.mat$'], 0);
save(analysisFile{1}, 'M', '-v7.3');
