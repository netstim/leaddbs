function analysisFile = ea_getGroupAnalysisFile(folder)
% Get group analysis file path from input dataset or group analysis folder.

folder = GetFullPath(folder);

if contains(folder, [filesep, 'derivatives', filesep, 'leadgroup', filesep]) % Input is group analysis folder
    dataset = regexp(folder, ['[^\' filesep ']+(?=\' filesep 'derivatives)'], 'match', 'once');
    pattern = ['^dataset-' dataset '_analysis-[^\W_]+\.mat$'];
    analysisFile = ea_regexpdir(folder, pattern, 0);
elseif isfolder(fullfile(folder, 'derivatives')) % Input is dataset root folder
    [~, dataset] = fileparts(erase(folder, filesep + textBoundary("end")));
    pattern = ['^dataset-' dataset '_analysis-[^\W_]+\.mat$'];
    analysisFile = ea_regexpdir(fullfile(folder, 'derivatives', 'leadgroup'), pattern, 1);
else
    error('Please specify either the dataset root folder or the group analysis folder!');
end

if ~isempty(analysisFile)
    % Choose one if multiple found
    if length(analysisFile) > 1
        [~, guid] = fileparts(fileparts(analysisFile));
        index = listdlg('PromptString', 'Select Group Analysis', 'ListString', guid, 'SelectionMode', 'single', 'CancelString', 'New Analysis');
        if isempty(index)
            analysisFile = [];
            return;
        else
            analysisFile = analysisFile(index);
        end
    end
    
    analysisFile = analysisFile{1};
else
    analysisFile = '';
end
