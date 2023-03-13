function ea_rename_pts(BIDSRoot, oldSubjId, newSubjId)
% Function to rename subj in a BIDS dataset

if ~iscell(oldSubjId)
    oldSubjId = {oldSubjId};
end

if ~iscell(newSubjId)
    newSubjId = {newSubjId};
end

if numel(oldSubjId) ~= numel(newSubjId)
    error('Length of old subjID doesn''t match length of new subjID!');
end

bids = BIDSFetcher(BIDSRoot);

if ~all(ismember(oldSubjId, bids.subjId))
    error('Not all the old subjID exist in the BIDS dataset!');
end

if numel(oldSubjId) ~= numel(unique(oldSubjId)) || numel(newSubjId) ~= numel(unique(newSubjId))
    error('Duplicated subjID found!');
end

if any(ismember(newSubjId, bids.subjId))
    error('New subjID conflicts with existing subjID in the BIDS dataset!');
end

for i = 1:numel(oldSubjId)
    old = oldSubjId{i};
    new = newSubjId{i};

    if isfolder(fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', old]))
        movefile(fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', old]), ...
                 fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', new]));
    end

    if isfolder(fullfile(BIDSRoot, 'rawdata', ['sub-', old]))
        movefile(fullfile(BIDSRoot, 'rawdata', ['sub-', old]), ...
                 fullfile(BIDSRoot, 'rawdata', ['sub-', new]));
    end

    if isfolder(fullfile(BIDSRoot, 'sourcedata', ['sub-', old]))
        movefile(fullfile(BIDSRoot, 'sourcedata', ['sub-', old]), ...
                 fullfile(BIDSRoot, 'sourcedata', ['sub-', new]));
    end

    oldFiles = ea_regexpdir(BIDSRoot, ['^sub-', old, '_']);
    newFiles = strrep(oldFiles, [filesep, 'sub-', old, '_'], [filesep, 'sub-', new, '_']);
    cellfun(@(src, dst) movefile(src, dst), oldFiles, newFiles);

    rawImageJson = fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', new], 'prefs', ['sub-', new, '_desc-rawimages.json']);
    if isfile(rawImageJson)
        json = fread(fopen(rawImageJson, 'rt'));
        json = strrep(char(json'), ['sub-', old, '_'], ['sub-', new, '_']);
        fwrite(fopen(rawImageJson, 'wt'), json);
    end

    statsFile = fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', new], ['sub-', new, '_desc-stats.mat']);
    if isfile(statsFile)
        load(statsFile, 'ea_stats');
        ea_stats.patname = strrep(ea_stats.patname, ['sub-', old], ['sub-', new]);
        save(statsFile, 'ea_stats');
    end

    statsBackupFile = fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', new], ['sub-', new, '_desc-statsbackup.mat']);
    if isfile(statsBackupFile)
        load(statsBackupFile, 'ea_stats');
        ea_stats.patname = strrep(ea_stats.patname, ['sub-', old], ['sub-', new]);
        save(statsBackupFile, 'ea_stats');
    end
end
