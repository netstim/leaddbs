function ea_generate_min_dataset(M, min_dataset_dir)

% Generates minimal datasets that contain electrode reconstructions
% and optionally VATs (including flipped)
% Such datasets allow VAT calculations in MNI space and
% can be processed in fiber filtering, neworkmapping and sweetspot tools 

% M               - lead-group file
% min_dataset_dir - directory to store the dataset, by default stored in LG 

% first we create varios directories to follow BIDS structure
if ~exist('min_dataset_dir', 'var')
    min_dataset_dir = [M.root,'Miniset'];
end

if exist(min_dataset_dir, 'dir')
    disp('Overwriting the old miniset')
    rmdir(min_dataset_dir, 's');
end

mkdir(min_dataset_dir);
mkdir([min_dataset_dir, '/derivatives']);
mkdir([min_dataset_dir, '/derivatives/leaddbs']);
mkdir([min_dataset_dir, '/derivatives/leadgroup']);
mkdir([min_dataset_dir, '/derivatives/leadgroup/', M.guid]);

% add the miniset flag
fid = fopen([min_dataset_dir, '/derivatives/leaddbs/Miniset_flag.json'],'w');
fid = fclose(fid);
%edit([min_dataset_dir, '/derivatives/leaddbs/Miniset_flag.json'])

% the leadgroup will be also stored, but detached
M_mini = M;
M_mini.root = [min_dataset_dir,'/derivatives/leadgroup/', M.guid];

for i = 1 : size(M.patient.list, 1)

    [filepath,patient_tag,ext] = fileparts(M.patient.list{i});
    mkdir([min_dataset_dir, '/derivatives/leaddbs/', patient_tag]);

    % here min_dataset_dir should be substituted with 'local directory/Miniset'
    M_mini.patient.list{i} = [min_dataset_dir, '/derivatives/leaddbs/', patient_tag];
    
    newPrefsFolder = [min_dataset_dir, '/derivatives/leaddbs/', patient_tag, '/prefs'];
    mkdir(newPrefsFolder);
    newReconstFolder = [min_dataset_dir, '/derivatives/leaddbs/', patient_tag, '/reconstruction'];
    mkdir(newReconstFolder);
    newStimFolder = [min_dataset_dir, '/derivatives/leaddbs/', patient_tag, '/stimulations/MNI152NLin2009bAsym/',M.S(i).label];
    mkdir(newStimFolder);
    
    % now we copy some essential files
    copyfile([M.patient.list{i},'/prefs/',patient_tag,'_desc-rawimages.json'], newPrefsFolder);
    copyfile([M.patient.list{i},'/reconstruction/',patient_tag,'_desc-reconstruction.mat'], newReconstFolder);  
    if isfile([M.patient.list{i},'/',patient_tag,'_desc-stats.mat'])
        copyfile([M.patient.list{i},'/',patient_tag,'_desc-stats.mat'], M_mini.patient.list{i});
    end

    % check if VATs were already computed 
    myVATs = dir(fullfile([M.patient.list{i}, '/stimulations/MNI152NLin2009bAsym/',M.S(i).label],[patient_tag,'_sim-efield_model*'])); %gets all mat files in struct
    % do not copy flipped VATs, they might be outdated!
    myVATs = myVATs(~endsWith({myVATs.name}, 'lippedFromLeft.nii'));
    myVATs = myVATs(~endsWith({myVATs.name}, 'lippedFromRight.nii'));

    if isempty(myVATs)
        ea_cprintf('CmdWinWarnings','No VATs were found for patient %s! Calculate them in Lead-Group', patient_tag);
    else
        for k = 1:length(myVATs)
            VAT_to_copy = fullfile(myVATs(k).folder, myVATs(k).name);
            copyfile(VAT_to_copy, newStimFolder);
        end
        
        % flip VTAs
        %ea_genflippedjointnii([newStimFolder,'/','*_hemi-R.nii'], [newStimFolder,'/','*_hemi-L.nii']); 

        if length(myVATs) == 2 % flipped will be overwritten
            % this order because L comes before R
            ea_genflippedjointnii([newStimFolder,'/',myVATs(2).name], [newStimFolder,'/',myVATs(1).name]);     
        elseif length(myVATs) == 1 && contains(myVATs(k).name, 'hemi-R')
            ea_genflippedjointnii([newStimFolder,'/',myVATs(2).name], 'no_VAT'); 
        elseif length(myVATs) == 1 && contains(myVATs(k).name, 'hemi-L')
            ea_genflippedjointnii('no_VAT', [newStimFolder,'/',myVATs(1).name]);
        else
            waitfor(msgbox('Wrong files were copied?'));
        end       
        
    end 

    % check if clinical scores are available
end

% save the updated lead-group file
M_orig = M;
M = M_mini;
[filepath,miniset_name,ext] = fileparts(min_dataset_dir);
save([M_mini.root,'/dataset-',miniset_name,'_analysis-',M_mini.guid,'.mat'], 'M')
M = M_orig;

disp('Done')

end
