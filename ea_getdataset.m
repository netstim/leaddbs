function uipatdir = ea_getdataset(options,handles)

try
    load([ea_prefsdir, filesep, 'ea_recentdatasets.mat'], 'recentfolders');
    startPath = fileparts(recentfolders{1});
catch
    % Use default location
    defaultLocation = fullfile(ea_gethome, 'Documents', 'LeadDBSDataset');
    if ~isfolder(defaultLocation)
        defaultLocationExisted = 0;
        ea_mkdir(defaultLocation);
    else
        defaultLocationExisted = 1;
    end
    startPath = defaultLocation;
end

uipatdir = uigetdir(startPath, 'Please choose dataset folder...');

if ~uipatdir
    return
end

% Delete created default dataset folder in case of choosing another folder
if exist('defaultLocation', 'var') && ~defaultLocationExisted && ~strcmp(uipatdir, defaultLocation)
    ea_delete(defaultLocation);
end

% Just in case "derivatives" or "derivatives/leaddbs" is chosen as dataset
% folder. This should not happen in the ideal case.
if endsWith(uipatdir, filesep + "derivatives") && ~isempty(ea_regexpdir(uipatdir, '^leaddbs$', 0, 'd'))
    uipatdir = erase(uipatdir, filesep + "derivatives" + lineBoundary("end"));
elseif endsWith(uipatdir, filesep + "derivatives" + filesep + "leaddbs")
    uipatdir = erase(uipatdir, filesep + "derivatives" + filesep + "leaddbs" + lineBoundary("end"));
end

[~, datasetName] = fileparts(uipatdir);
if ~isempty(regexp(datasetName, '[\W_]', 'once'))
    ea_error(sprintf('Should only contain alphanumeric characters but "%s" provided!\n', datasetName), title = 'Please adapt dataset folder name', simpleStack = 1);
end

ea_mkdir(fullfile(uipatdir, 'derivatives', 'leaddbs'));
ea_mkdir(fullfile(uipatdir, 'rawdata'));
ea_mkdir(fullfile(uipatdir, 'sourcedata'));

if ~isfile(fullfile(uipatdir, 'dataset_description.json'))
    ea_cprintf('CmdWinWarnings', 'Could not find dataset description file, generating one now...\n');
    ea_generate_datasetDescription(uipatdir, 'root_folder');
end

if exist('handles','var')
    ea_load_pts(handles, uipatdir);

    if isfield(handles, 'atlassetpopup') % not present in connectome mapper
        atlasset=get(handles.atlassetpopup, 'String');
        atlasset=atlasset{get(handles.atlassetpopup, 'Value')};

        ea_listatlassets(options, handles, get(handles.vizspacepopup, 'Value'), atlasset);
    end
end
