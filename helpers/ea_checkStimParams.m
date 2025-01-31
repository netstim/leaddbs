function output = ea_checkStimParams(stimFile)
% Check if "S" is compatible with the latest format. Fix when necessary.
%
% Input can be stimparameters file or group analysis file.
% Output will be "S" or "M" (with "M.S" updated).

updated = 0;

stimFile = GetFullPath(stimFile);
vars = who('-file', stimFile);
if ismember('M', vars)
    load(stimFile, 'M');
    S = M.S;
elseif ismember('S', vars)
    load(stimFile, 'S');
end

% Return and do not patch fields when initialzing S
if isempty(S)
    if exist('M', 'var')
        output = M;
    else
        output = S;
    end
    return;
end

if ~isfield(S, 'sources')
    [S.sources] = deal(1:4);
    updated = 1;
end

if ~isfield(S, 'volume')
    [S.volume] = deal([]);
    updated = 1;
end

if ~isfield(S, 'numContacts')
    if exist('M', 'var')
        numContacts = cell(1, length(M.elstruct));
        for i=1:length(M.elstruct)
            options.elmodel = M.elstruct(i).elmodel;
            options = ea_resolve_elspec(options);
            numContacts{i} = options.elspec.numContacts;
        end
        [S.numContacts] = deal(numContacts{:});
    else
        subjFolder = regexp(stimFile, ['.+derivatives\' filesep 'leaddbs\' filesep 'sub-[^\W_]+'], 'match');
        uiprefFile = ea_regexpdir(fullfile(subjFolder{1}, 'prefs'), '_desc-uiprefs\.mat', 0, 'f');
        if ~isempty(uiprefFile) && isfile(uiprefFile{1}) % Check uiprefs first
            options = ea_resolve_elspec(load(uiprefFile{1}, 'elmodel'));
        else
            reconFile = ea_regexpdir(fullfile(subjFolder{1}, 'reconstruction'), '_desc-reconstruction\.mat', 0, 'f');
            if ~isempty(reconFile) && isfile(reconFile{1}) % Check reconstruction
                load(reconFile{1}, 'reco');
                options.elmodel = ea_get_first_notempty_elmodel(reco.props);
                options = ea_resolve_elspec(options);
            else
                ea_error('uiprefs and reconstruction are both not present! Failed to detect electrode.', showdlg=0, simpleStack=1);
            end
        end
        S.numContacts = options.elspec.numContacts;
    end

    updated = 1;
end

if ~isfield(S, 'ver')
    % Further fix contact field
    S = arrayfun(@adaptContactField, S);

    [S.ver] = deal('2.0');
    updated = 1;
end

if exist('M', 'var')
    M.S = S;
    output = M;
else
    output = S;
end

if updated
    if exist('M', 'var')
        save(stimFile, 'M');
    else
        save(stimFile, 'S');
    end
end


function newS = adaptContactField(S)
%% Initialize new S first
newS = S;
% Right sources
for source=1:4
    newS.(['Rs',num2str(source)]) = struct;
    for k=1:S.numContacts
        newS.(['Rs',num2str(source)]).(['k',num2str(k)]).perc=0;
        newS.(['Rs',num2str(source)]).(['k',num2str(k)]).pol=0;
        newS.(['Rs',num2str(source)]).(['k',num2str(k)]).imp=1;
    end
    newS.(['Rs',num2str(source)]).case.perc = 100;
    newS.(['Rs',num2str(source)]).case.pol = 2;
    newS.(['Rs',num2str(source)]).amp = 0;
    newS.(['Rs',num2str(source)]).va = 2;
    newS.(['Rs',num2str(source)]).pulseWidth = 60;
end

% Left sources
for source=1:4
    newS.(['Ls',num2str(source)]) = struct;
    for k=1:S.numContacts
        newS.(['Ls',num2str(source)]).(['k',num2str(k)]).perc=0;
        newS.(['Ls',num2str(source)]).(['k',num2str(k)]).pol=0;
        newS.(['Ls',num2str(source)]).(['k',num2str(k)]).imp=1;
    end
    newS.(['Ls',num2str(source)]).case.perc = 100;
    newS.(['Ls',num2str(source)]).case.pol = 2;
    newS.(['Ls',num2str(source)]).amp = 0;
    newS.(['Ls',num2str(source)]).va = 2;
    newS.(['Ls',num2str(source)]).pulseWidth = 60;
end

%% Copy old stimulations
% Right sources
for source=1:4
    for k=1:min(S.numContacts,8)
        newS.(['Rs',num2str(source)]).(['k',num2str(k)]) = S.(['Rs',num2str(source)]).(['k',num2str(k-1)]);
    end
    newS.(['Rs',num2str(source)]).case = S.(['Rs',num2str(source)]).case;
    newS.(['Rs',num2str(source)]).amp = S.(['Rs',num2str(source)]).amp;
    newS.(['Rs',num2str(source)]).va = S.(['Rs',num2str(source)]).va;
    if isfield(S.(['Rs',num2str(source)]), 'pulseWidth')
        newS.(['Rs',num2str(source)]).pulseWidth = S.(['Rs',num2str(source)]).pulseWidth;
    end
end

% Left sources
for source=1:4
    for k=1:min(S.numContacts,8)
        newS.(['Ls',num2str(source)]).(['k',num2str(k)]) = S.(['Ls',num2str(source)]).(['k',num2str(k+7)]);
    end
    newS.(['Ls',num2str(source)]).case = S.(['Ls',num2str(source)]).case;
    newS.(['Ls',num2str(source)]).amp = S.(['Ls',num2str(source)]).amp;
    newS.(['Ls',num2str(source)]).va = S.(['Ls',num2str(source)]).va;
    if isfield(S.(['Ls',num2str(source)]), 'pulseWidth')
        newS.(['Ls',num2str(source)]).pulseWidth = S.(['Ls',num2str(source)]).pulseWidth;
    end
end
