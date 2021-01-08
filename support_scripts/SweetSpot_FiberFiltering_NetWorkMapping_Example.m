clear
lead path % add Lead-DBS to the path

resultfig = ea_mnifigure; % Create empty 3D viewer figure

M.pseudoM = 1; % Declare this is a pseudo-M struct, i.e. not a real lead group file
M.ROI.list={'/path/to/first_nifti.nii' % enter paths to ROI here
'/path/to/second_nifti.nii'
'/path/to/third_nifti.nii'
'/path/to/fourth_nifti.nii'
'/path/to/fifth_nifti.nii'
'/path/to/sixth_nifti.nii'
'/path/to/seven_nifti.nii'
};

M.ROI.group=ones(length(M.ROI.list),1);

M.clinical.labels={'Improvement','Covariate_1','Covariate_2'}; % how will variables be called
M.clinical.vars{1}=[1
    2
    3
    4
    5
    6
    7]; % enter a variable of interest - entries correspond to nifti files
M.clinical.vars{2}=[1.4
    2.5
    0.3
    4.4
    5.5
    0.6
    8.7]; % enter as many variables of interest as you want
M.clinical.vars{3}=[1.3
    2.5
    0.3
    10.4
    1.5
    2.6
    9.7]; % enter as many variables of interest as you want

M.guid='My_Analysis'; % give your analysis a name

save('Analysis_Input_Data.mat','M'); % store data of analysis to file

% Open up the Sweetspot Explorer
ea_sweetspotexplorer(fullfile(pwd,'Analysis_Input_Data.mat'),resultfig);

% Open up the Network Mapping Explorer
ea_networkmappingexplorer(fullfile(pwd,'Analysis_Input_Data.mat'),resultfig);

% Open up the Fiber Filtering Explorer
ea_discfiberexplorer(fullfile(pwd,'Analysis_Input_Data.mat'),resultfig);





