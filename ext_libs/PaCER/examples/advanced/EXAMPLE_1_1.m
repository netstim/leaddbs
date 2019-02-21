%% PaCER/EXAMPLE_1_1.m - Example file for using PaCER - Precise and Convenient Electrode Reconstruction
%                       to reconstruct precise curved electrode trajectories
%                       and contact locations from post operative CT imaging
%                       of DBS electrodes
%
% Andreas Husch
% Centre Hospitalier de Luxembourg (CHL) / Luxembourg Centre for Systems
% Biomedicine (LCSB), University of Luxembourg
% mail@andreashusch.de / husch.andreas@chl.lu
% 2017
%
%% Make sure all PaCER files (including subdirectories) are in your MATLAB path! 
SETUP_PACER

%% Define a post OP CT nifti file & load it into matlab for further use
POSTOPCT_FILENAME = 'native_nifti_data/DEMO_postopCT.native.nii.gz'; 
niiCT = NiftiMod(POSTOPCT_FILENAME);

%% Run PACER on the post OP CT
% the result "elecModels" is a set of electrodes objects reconstructed from the CT, i.e. elecModels{1} is the first electrode found, elecModels{2} the second
elecModels = PaCER(niiCT); 

%% Plot the electrodes
figure;
elecModels{1}.initPlot3D(gca);
elecModels{2}.initPlot3D(gca);

%% Continuing from EXAMPLE_1 here:

%%  Add som some MPR (axial, coronal, sagital) of the CT
%  you can scroll the MPR by clicking a slice and turning the mouse wheel!
mpr = createSimpleMPRWorldCoordinates(niiCT);
colormap gray;

%% ****** All subsequent code customizes the plot a bit for nicer display. Note that PaCER ******
% electrode models feature automatic plot updates when changing key parameters on the fly 
% (like changinge electrode color, transparency or contact display mode)

%% *** Some custom PaCER toolbox plot features ***
%% change the CT colormap using a non-monotonic tone-mapping
CTColormap();

%% *** Some plot features of the PolynomialElectrodeModel class (automatic plot updates) ***
%% Change primary color of the first electrode
elecModels{1}.ELECTRODE_COLOR = rgb('LimeGreen');

%% Make the second electrodes more transparent (default ALPHA=0.9)
elecModels{2}.ALPHA = 0.5;

%% Display contact locations as expected from the electrode type/model instead as detected
%  (this might by useful if either your slice thickness is to coarse for
%  precise single contact detection OR if you have a good single contact 
%  detection and want to compare to the expected locations from the model)
elecModels{1}.useDetectedContactPositions = false;

%% *** Some standard MATLAB plot option ***
axis tight
xlabel('X')
ylabel('Y')
zlabel('Z')
lighting gouraud
camlight headlight;
camlight left;
grid on
view(-2,2)
camtarget(elecModels{1}.getEstTipPos);
cameratoolbar
cameratoolbar('nomode')
camzoom(2)
% use the "orbit camera tool" (upper left corner) to rotate view.
% deactivate the camera to move an mpr slice by clicking it and turning the
% mouse wheel
