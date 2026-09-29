function pairs = ea_lc_checkreg_pairs(options, checkDMRI, checkFMRI)
% Resolve BIDS/legacy Lead-Connectome registration image pairs.

arguments
    options struct
    checkDMRI (1,1) logical = false
    checkFMRI (1,1) logical = false
end

directory = fullfile(options.root, options.patientname);
pairs = struct('moving', {}, 'fixed', {}, 'tag', {});

if checkDMRI && isfile(fullfile(directory, 'coregistration', 'log', 'ea_lc_b0_coreg.mat'))
    record = ea_lc_coreg_record(options);
    pairs(end+1) = makePair(record.output, record.anchor, 'Lead-Connectome dMRI: diffusion & anatomy');
    checkDMRI = false;
end

if checkDMRI
    [options, ~, b0] = ea_lc_resolve_b0(options);
    anat = resolveFile(directory, options.prefs, 'prenii_unnormalized');
    moving = '';
    if ~isempty(b0) && ~isempty(anat)
        [~, b0Name] = ea_niifileparts(b0);
        [~, anatName] = ea_niifileparts(anat);
        candidates = {
            fullfile(directory, 'preprocessing', 'anat', [b0Name, '2', anatName, '.nii'])
            fullfile(directory, 'coregistration', 'anat', ...
                [options.patientname, '_space-anchorNative_dwi_fa.nii'])
        };
        moving = firstExisting(candidates);
        if isempty(moving)
            moving = findMatchingNifti(fullfile(directory, 'coregistration', 'dwi'), ...
                {b0Name, anatName});
        end
    end
    if ~isempty(moving)
        % B0-to-anatomy products are in anatomical space. Older tracking-mask
        % products are anatomy-to-B0 and therefore match B0 space instead.
        if contains(moving, [filesep, 'coregistration', filesep, 'dwi', filesep])
            fixed = b0;
        else
            fixed = anat;
        end
        pairs(end+1) = makePair(moving, fixed, 'Lead-Connectome dMRI: diffusion & anatomy');
    end
end

if checkFMRI
    anat = resolveFile(directory, options.prefs, 'prenii_unnormalized');
    restFiles = findRestFiles(directory, options.prefs);
    if isempty(anat)
        restFiles = {};
        anatName = '';
    else
        [~, anatName] = ea_niifileparts(anat);
    end
    for i = 1:numel(restFiles)
        rest = restFiles{i};
        [restStem, restName] = ea_niifileparts(rest);
        restDir = fileparts(restStem);
        meanRest = fullfile(restDir, ['mean', restName, '.nii']);
        registeredAnat = fullfile(restDir, ['r', restName, '_', anatName, '.nii']);
        if isfile(meanRest) && isfile(registeredAnat)
            pairs(end+1) = makePair(registeredAnat, meanRest, ...
                ['Lead-Connectome fMRI: ', restName, ' & anatomy']); %#ok<AGROW>
        end
    end
end


function path = resolveFile(directory, prefs, field)
path = '';
if isfield(prefs, field) && ~isempty(prefs.(field))
    candidate = prefs.(field);
    if ~isfile(candidate)
        candidate = fullfile(directory, candidate);
    end
    if isfile(candidate)
        path = candidate;
    end
end


function path = firstExisting(candidates)
path = '';
for i = 1:numel(candidates)
    if isfile(candidates{i})
        path = candidates{i};
        return;
    end
end


function path = findMatchingNifti(folder, tokens)
path = '';
if ~isfolder(folder)
    return;
end
files = [dir(fullfile(folder, '*.nii')); dir(fullfile(folder, '*.nii.gz'))];
for i = 1:numel(files)
    if all(cellfun(@(token) contains(files(i).name, token), tokens))
        path = fullfile(files(i).folder, files(i).name);
        return;
    end
end


function files = findRestFiles(directory, prefs)
files = {};
if isfield(prefs, 'rest') && ~isempty(prefs.rest)
    rests = cellstr(prefs.rest);
    for i = 1:numel(rests)
        candidate = rests{i};
        if ~isfile(candidate)
            candidate = fullfile(directory, candidate);
        end
        if isfile(candidate)
            files{end+1} = candidate; %#ok<AGROW>
        end
    end
end
if isempty(files)
    found = [dir(fullfile(directory, 'preprocessing', 'func', '*_bold.nii')); ...
             dir(fullfile(directory, 'preprocessing', 'func', '*_bold.nii.gz'))];
    keep = ~cellfun(@(name) startsWith(name, {'r', 'sr', 'mean', 'hdmean'}), ...
        {found.name});
    found = found(keep);
    files = arrayfun(@(f) fullfile(f.folder, f.name), found, 'UniformOutput', false);
end


function pair = makePair(moving, fixed, tag)
pair = struct('moving', moving, 'fixed', fixed, 'tag', tag);
