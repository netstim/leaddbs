function ea_generate_min_dataset(M, min_dataset_dir)
% Generates minimal datasets that contain electrode reconstructions
% and optionally VATs (including flipped)
% Such datasets allow VAT calculations in MNI space and
% can be processed in fiber filtering, neworkmapping and sweetspot tools 

% M               - lead-group file
% min_dataset_dir - directory to store the dataset, by default stored in LG 

% first we create varios directories to follow BIDS structure
if ~exist('min_dataset_dir', 'var')
    min_dataset_dir = fullfile(M.root, 'Miniset');
end

if exist(min_dataset_dir, 'dir')
    disp('Overwriting the old miniset!')
    rmdir(min_dataset_dir, 's');
end

mkdir(fullfile(min_dataset_dir, 'derivatives', 'leaddbs'));
mkdir(fullfile(min_dataset_dir, 'derivatives', 'leadgroup', M.guid));

% add the miniset flag
[~, minisetName] = fileparts(min_dataset_dir);
miniset.name = minisetName;
miniset.numSubj = length(M.patient.list);
miniset.recon = 1;
miniset.vta = 0;
miniset.stats = 0;
miniset.groupAnalysis = M.guid;
miniset.timeStamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

% the leadgroup will be also stored, but detached
M_mini = M;
M_mini.root = fullfile(min_dataset_dir, 'derivatives', 'leadgroup', M.guid, filesep);

for i = 1 : size(M.patient.list, 1)
    [~, patient_tag] = fileparts(M.patient.list{i});
    fprintf('Subject %d/%d: %s ...\n', i, size(M.patient.list, 1), patient_tag);

    newPatientFolder = fullfile(min_dataset_dir, 'derivatives', 'leaddbs', patient_tag);
    ea_mkdir(newPatientFolder);

    % here min_dataset_dir should be substituted with 'local directory/Miniset'
    M_mini.patient.list{i} = newPatientFolder;
    
    newPrefsFolder = fullfile(newPatientFolder, 'prefs');
    ea_mkdir(newPrefsFolder);
    newReconstFolder = fullfile(newPatientFolder, 'reconstruction');
    ea_mkdir(newReconstFolder);
    newStimFolder = fullfile(newPatientFolder, 'stimulations', ea_getspace, M.S(i).label);
    ea_mkdir(newStimFolder);
    
    % now we copy some essential files
    if isfile(fullfile(M.patient.list{i}, 'prefs', [patient_tag, '_desc-rawimages.json']))
        ea_cprintf('CmdWinWarnings', 'No rawimages.json found for patient %s!\n', patient_tag);
        copyfile(fullfile(M.patient.list{i}, 'prefs', [patient_tag, '_desc-rawimages.json']), newPrefsFolder);
    end
    %copyfile(fullfile(M.patient.list{i}, 'prefs', [patient_tag, '_desc-uiprefs.mat']), newPrefsFolder);
    copyfile(fullfile(M.patient.list{i}, 'reconstruction', [patient_tag, '_desc-reconstruction.mat']), newReconstFolder);  
    if isfile(fullfile(M.patient.list{i}, [patient_tag, '_desc-stats.mat']))
        miniset.stats = 1;
        copyfile(fullfile(M.patient.list{i}, [patient_tag, '_desc-stats.mat']), M_mini.patient.list{i});
    end

    % Copy stimulation parameters
    stimParams = fullfile(M.patient.list{i}, 'stimulations', ea_getspace, M.S(i).label, [patient_tag,'_desc-stimparameters.mat']);
    if isfile(stimParams)
        copyfile(stimParams, newStimFolder);
    end

    % check if VATs were already computed 
    myVATs = dir(fullfile(M.patient.list{i}, 'stimulations', ea_getspace, M.S(i).label, [patient_tag,'_sim-*'])); %gets all mat files in struct
    % do not copy flipped VATs, they might be outdated!
    myVATs = myVATs(~contains({myVATs.name}, 'flippedFrom', 'IgnoreCase', true));
    % Ignore mapper results
    myVATs = myVATs(~contains({myVATs.name}, '_seed-'));

    if isempty(myVATs)
        ea_cprintf('CmdWinWarnings','No VATs were found for patient %s! Calculate them in Lead-Group \n', patient_tag);
    else
        miniset.vta = 1;
        for k = 1:length(myVATs)
            VAT_to_copy = fullfile(myVATs(k).folder, myVATs(k).name);
            copyfile(VAT_to_copy, newStimFolder);
        end       
    end 

    % check if clinical scores are available
end

% save the updated lead-group file
M = M_mini;
save(fullfile(M_mini.root, ['dataset-', minisetName, '_analysis-', M_mini.guid, '.mat']), 'M')

% Save miniset json
savejson('', miniset, fullfile(min_dataset_dir, 'miniset.json'));

disp('Done')
