%% PaCER/EXAMPLE_2_1.m - Longitudinal Datasets and Plot of Segmentations and simplified VTA 
%
% Andreas Husch
% Centre Hospitalier de Luxembourg (CHL) / Luxembourg Centre for Systems
% Biomedicine (LCSB), University of Luxembourg
% mail@andreashusch.de / husch.andreas@chl.lu
% 2017
%
%% Make sure all PaCER files (including subdirectories) are in your MATLAB path!
SETUP_PACER % SETUP_PACER has to be run *once* to setup the PaCER path, further runs are not required

%% Load two co-registred post OP CT files aquired at different time points
POSTOPCT_FILENAME_1 = 'pre-registered_nifti_data/DEMO_postopCT.preopCTspace.nii.gz'; 
POSTOPCT_FILENAME_2 = 'pre-registered_nifti_data/DEMO_followupCT.preopCTspace.nii.gz'; 
niiPostOP2 = NiftiMod(POSTOPCT_FILENAME_2); % load the second file

PREOPT1_FILENAME1 = 'pre-registered_nifti_data/DEMO_T1.preopCTspace.nii.gz'; % registered to preop
niiT1 = NiftiMod(PREOPT1_FILENAME1);

%% Run PACER
elecModels1 = PaCER(POSTOPCT_FILENAME_1); % PaCER can work with a file NAME as input...
elecModels2 = PaCER(POSTOPCT_FILENAME_2); % ...or with a NiftiMod object

%% Plot MPR of preop T1 MRI
figure('Name', 'Example 2.1');
mpr = createSimpleMPRWorldCoordinates(niiT1); %scroll MPR planes by clicking a plane and turing the mouse wheel
colormap gray
hold on;

%% Plot electrodes
elecModels1{1}.ELECTRODE_COLOR = rgb('Crimson');
elecModels1{2}.ELECTRODE_COLOR = rgb('Crimson');
elecModels1{1}.plot3D;
elecModels1{2}.plot3D;

elecModels2{1}.plot3D;
elecModels2{2}.plot3D;


%% Load Manual STN Segmentations
stnL = NiftiSeg('pre-registered_nifti_data/STN_lh.preopCTspace.nii.gz');
stnR = NiftiSeg('pre-registered_nifti_data/STN_rh.preopCTspace.nii.gz');

%% Plot the Segemntations
stnL.plot3D;
stnR.plot3D;

%% Change Colors 
stnL.setLabelColorByIdx(2, rgb('LimeGreen')); %change by index. If a label file is given changes can also be made by names (see EXAMPLE_4)
stnR.setLabelColorByIdx(2, rgb('LimeGreen'));

%% Add and display a simplied monopolar field simulation for electrodeModel1{1} (Mädler / Coenen Model)
%  (added to contact 1 [note: we use matlab counting terminolgy, i.e. starting at *1* at
%  the tip, might reflect "contact 0" in some vendor terminolgy])
fs = SimpleFieldSim3D(elecModels1{1});
fs.plot3D

%% Switch field model to contact 2 on the fly, try others..
elecModels1{1}.activeContact = 2;

%% Increase voltage to 2.0V, try other values..
fs.voltage = 2.0;

%% %Display a GUI to tune field simulation paratmeters graphically on the fly 
fs.getConfigWindow()
