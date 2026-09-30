function modality = ea_getmodality(BIDSFilePath)
% Extract image modality from BIDS file path

if ~iscell(BIDSFilePath)
    wasChar = 1;
    BIDSFilePath = {BIDSFilePath};
else
    wasChar = 0;
end

modality = cell(size(BIDSFilePath));

for i=1:length(BIDSFilePath)
    try
        parsedStruct = parseBIDSFilePath(BIDSFilePath{i});
        hasAcq = isfield(parsedStruct, 'acq')    && ~isempty(parsedStruct.acq);
        hasSuf = isfield(parsedStruct, 'suffix') && ~isempty(parsedStruct.suffix);
        isCT   = hasSuf && strcmp(parsedStruct.suffix, 'CT');

        if hasAcq && hasSuf && ~isCT
            modality{i} = [parsedStruct.acq '_' parsedStruct.suffix];
        elseif hasSuf
            modality{i} = parsedStruct.suffix;
        elseif hasAcq
            modality{i} = parsedStruct.acq;
        end
    catch
        % Fallback for filenames that don't strictly conform to BIDS:
        % strip all key-value entities (e.g. 'sub-XX_', 'ses-preop_',
        % 'desc-preproc_') and use whatever remains as the modality token.
        [~, fname] = fileparts(BIDSFilePath{i});
        modality{i} = regexprep(fname, '[a-zA-Z]+-[^\W_]+_', '');
    end
end

if wasChar
    modality = modality{1};
end
