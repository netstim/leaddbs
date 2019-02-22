%% PaCER/EXAMPLE_4.m - Demonstrate display atlas data with PaCER reconstructions from CT's previouly transformed to atlas space 
%
% Andreas Husch
% Centre Hospitalier de Luxembourg (CHL) / Luxembourg Centre for Systems
% Biomedicine (LCSB), University of Luxembourg
% mail@andreashusch.de / husch.andreas@chl.lu
% 2017
%
%% Make sure all PaCER files (including subdirectories) are in your MATLAB path!
SETUP_PACER % SETUP_PACER has to be run *once* to setup the PaCER path, further runs are not required

%% Load a post OP CT files in deep7T / MNI Space
POSTOPCT_FILENAME = 'atlas_coreg_data/DEMO_postopCT.template-space.nii.gz'; % transformed to deep7T atlas space
PREOPT1_FILENAME = 'atlas_coreg_data/DEMO_T1.template-space.nii.gz'; % registered to postop_ct1
atlas = NiftiSeg('aseg_deepBrainNuclei.nii.gz', 'labelFilePath', 'Deepbrain_7T.label');

niiT1 = NiftiMod(PREOPT1_FILENAME);
%% Run PACER
elecModels1 = PaCER(POSTOPCT_FILENAME, 'electrodeType', 'Medtronic 3389'); % we tell PaCER the electrode type manually, as the SNR is reduced due to the  transformation of the post OP CT to atlas space
%elecModels = PaCER(POSTOPCT_FILENAME, 'finalDegree', 1, 'electrodeType', 'Medtronic 3389'); % we tell PaCER the electrode type manually, as the SNR is reduced due to the  transformation of the post OP CT to atlas space

%% Plot MPR
figure('Name', 'Example 4');
mpr = createSimpleMPRWorldCoordinates(niiT1); %scroll MPR planes by clicking a plane and turing the mouse wheel
colormap gray
hold on;
camzoom(2);

%% Plot electrodes
elecModels{1}.plot3D;
elecModels{2}.plot3D;

%% Plot Atlas Structures
atlas.plot3D;
atlas.setLabelColorByName(atlas.presentLabelNames, rgb('Goldenrod'));
atlas.setLabelColorByName({'STN_L', 'STN_R'}, rgb('LimeGreen'));
