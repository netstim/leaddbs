function S = ea_checkStimParams(S)
% Check if "S" is compatible with the latest format. Fix when necessary.

updated = 0;

if isfile(S)
    stimFile = GetFullPath(S);
    var = who('-file', stimFile);
    if ismember('M', var)
        load(stimFile, 'M');
        S = M.S;
    elseif ismember('S', var)
        load(stimFile, 'S');
    end
end

if ~isfield(S, 'sources')
    [S.sources] = deal(1:4);
    updated = 1;
end

if ~isfield(S, 'volume')
    [S.volume] = deal([]);
    updated = 1;
end

if ~isfield(S, 'ver')
    [S.ver] = deal('2.0');
    updated = 1;
end

if ~isfield(S, 'numel')
    if exist('M', 'var')
        numel = cell(1, length(M.elstruct));
        for i=1:length(M.elstruct)
            options.elmodel = M.elstruct(i).elmodel;
            options = ea_resolve_elspec(options);
            numel{i} = options.elspec.numel;
        end
        [S.numel] = deal(numel{:});
    else
        subjFolder = regexp(stimFile, ['.+derivatives\' filesep 'leaddbs\' filesep 'sub-[^\W_]+'], 'match');
        uiprefFile = ea_regexpdir(fullfile(subjFolder{1}, 'prefs'), '_desc-uiprefs\.mat', 0, 'f');
        if isfile(uiprefFile{1}) % Check uiprefs first 
            options = ea_resolve_elspec(load(uiprefFile{1}, 'elmodel'));
        else
            reconFile = ea_regexpdir(fullfile(subjFolder{1}, 'reconstruction'), '_desc-reconstruction\.mat', 0, 'f');
            if isfile(reconFile{1}) % Check reconstruction
                load(reconFile{1}, 'reco');
                options.elmodel = ea_get_first_notempty_elmodel(reco.props);
                options = ea_resolve_elspec(options);
            else
                ea_error('uiprefs and reconstruction are both not present! Failed to detect electrode.', showdlg=0, simpleStack=1);
            end
        end
        S.numel = options.elspec.numel;
    end
    
    updated = 1;
end

if updated && exist('filePath', 'var')
    if exist('M', 'var')
        M.S = S;
        save(stimFile, 'M');
    else
        save(stimFile, 'S');
    end
end
