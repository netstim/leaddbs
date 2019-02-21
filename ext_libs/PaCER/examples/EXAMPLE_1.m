%% PaCER/EXAMPLE.m - Example file for using PaCER - Precise and Convenient Electrode Reconstruction
%                    to reconstruct precise curved electrode trajectories
%                    and contact locations from post operative CT imaging
%                    of DBS electrodes
%
% Andreas Husch
% Centre Hospitalier de Luxembourg (CHL) / Luxembourg Centre for Systems
% Biomedicine (LCSB), University of Luxembourg
% mail@andreashusch.de / husch.andreas@chl.lu
% 2017
%
%% Make sure all PaCER files (including subdirectories) are in your MATLAB path!
SETUP_PACER % SETUP_PACER has to be run *once* to setup the PaCER path, further runs are not required

%% Define a post OP CT nifti file & load it into matlab for further use
POSTOPCT_FILENAME = 'post_op_ct.nii.gz'; 
niiCT = NiftiMod(POSTOPCT_FILENAME);

%% Run PACER on the post OP CT
% the result "elecModels" is a set of electrodes objects reconstructed from the CT, 
% i.e. elecModels{1} corresponds to the first electrode found, elecModels{2} to the second
elecModels = PaCER(niiCT); % CT is the only required input, most convienient ;-)

%% Plot the electrodes
figure('Name', 'Example 1');
elecModels{1}.initPlot3D(gca);
elecModels{2}.initPlot3D(gca);

%% See EXAMPLE_1_1 for adding further visual elements for "standalone use" of PaCER
% alternativley, just use the results in elecModels with your favorite DBS
% toolbox :-)
