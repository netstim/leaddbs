function ea_import_patient(dataset, patientID, images, electrodeModel)
arguments
    dataset         {mustBeTextScalar}  % dataset folder
    patientID       {mustBeTextScalar}  % patient ID
    images          {mustBeText}        % patient images
    electrodeModel  {mustBeText} = ''   % electrode model
end

dataset = GetFullPath(dataset);

%% make dataset dirs
derivativeFolder = fullfile(dataset, 'derivatives', 'leaddbs');
rawDataFolder = fullfile(dataset, 'rawdata');
ea_mkdir(derivativeFolder);
ea_mkdir(rawDataFolder);
ea_mkdir(fullfile(dataset, 'sourcedata'));

if ~isfile(fullfile(dataset, 'dataset_description.json'))
    ea_generate_datasetDescription(dataset);
end

subPrefix = ['sub-', patientID];
patientDerivativeFolder = fullfile(derivativeFolder, subPrefix);
preopRawdataFolder = fullfile(rawDataFolder, subPrefix, 'ses-preop', 'anat');
postopRawdataFolder = fullfile(rawDataFolder, subPrefix, 'ses-postop', 'anat');
ea_mkdir(patientDerivativeFolder);
ea_mkdir(fullfile(patientDerivativeFolder, 'prefs'));
ea_mkdir(preopRawdataFolder);
ea_mkdir(postopRawdataFolder);

ea_cprintf('*Comments', 'Preparing patient folder %s ...\n', subPrefix);

%% copy raw data
if ischar(images)
    images = {images};
end

for i=1:numel(images)
    if contains(images{i}, {'T1', 'T2', 'FLAIR', 'FGATIR'}, 'IgnoreCase', 1)
        acqTag = ea_checkacq(images{i});
        if contains(images{i}, 'T1', 'IgnoreCase', 1)
            ea_cprintf('*Comments', 'Found pre-operative T1w ...\n');
            acqModality = [acqTag, '_T1w'];
        elseif contains(images{i}, 'T2', 'IgnoreCase', 1)
            ea_cprintf('*Comments', 'Found pre-operative T2w ...\n');
            acqModality = [acqTag, '_T2w'];
        elseif contains(images{i}, 'FLAIR', 'IgnoreCase', 1)
            ea_cprintf('*Comments', 'Found pre-operative FLAIR ...\n');
            acqModality = [acqTag, '_FLAIR'];
        elseif contains(images{i}, 'FGATIR', 'IgnoreCase', 1)
            ea_cprintf('*Comments', 'Found pre-operative FGATIR ...\n');
            acqModality = [acqTag, '_FGATIR'];
        end
        preopRawImage = fullfile(preopRawdataFolder, [subPrefix, '_ses-preop_acq-', acqModality, '.nii']);
        copyfile(images{i}, preopRawImage);
        gzip(preopRawImage);
        delete(preopRawImage);
    elseif contains(images{i}, 'CT', 'IgnoreCase', 1)
        ea_cprintf('*Comments', 'Found post-operative CT ...\n');
        postopRawImage = fullfile(postopRawdataFolder, [subPrefix, '_ses-postop_CT.nii']);
        copyfile(images{i}, postopRawImage);
        gzip(postopRawImage);
        delete(postopRawImage);
    end
end

%% setup rawimage.json
ea_genrawimagesjson(dataset, patientID);

%% Save uiprefs, set electrode model
uiprefsFile = fullfile(patientDerivativeFolder, 'prefs', [subPrefix, '_desc-uiprefs.mat']);
if isfile(uiprefsFile)
    uiprefs = load(uiprefsFile);
else
    uiprefs = loadjson(fullfile(ea_getearoot, 'common', 'uiprefs.json'));
    uiprefs.earoot = ea_getearoot;
end

if ~isempty(electrodeModel)
    uiprefs.elmodel = electrodeModel;
    uiprefs.elmodeln = find(strcmp(ea_resolve_elspec, electrodeModel));
end

save(uiprefsFile, '-struct', 'uiprefs');
