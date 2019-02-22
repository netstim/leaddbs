%% PaCER/EXAMPLE_5_1.m - Continuing Example 5 by adding fibertracts from MRtrix and projecting a VTA to them
%
% Andreas Husch
% Centre Hospitalier de Luxembourg (CHL) / Luxembourg Centre for Systems
% Biomedicine (LCSB), University of Luxembourg
% mail@andreashusch.de / husch.andreas@chl.lu
%
% Mikkel V. Petersen
% Center of Functionally Integrative Neuroscience (CFIN)
% Aarhus University Hospital
% mikkel.petersen@cfin.au.dk
%
% 2017
%
%% Make sure all PaCER files (including subdirectories) are in your MATLAB path! 
SETUP_PACER

%% Define post OP CT and T1 nifti files & load them into matlab for further use
POSTOPCT_FILENAME = 'native_nifti_data/DEMO_postopCT.native.nii.gz'; 
niiCT = NiftiMod(POSTOPCT_FILENAME);

T1_FILENAME = 'native_nifti_data/DEMO_T1.native.nii.gz';
niiT1 = NiftiMod(T1_FILENAME);

%% Define subcortical atlas nonlinearly transformed (FNIRT) to T1 space
ASEG_T1 = 'ASEG_nonlin.t1space.nii.gz';
atlas = NiftiSeg(ASEG_T1, 'labelFilePath', 'Deepbrain_7T.label');

%% Define and load FLIRT transformation matrix
FLIRT_TRANSFORM = 'native_nifti_data/flirt_transform.postopCT_to_T1.mat';
FlirtMat = load(FLIRT_TRANSFORM,'-ascii');

%% Run PACER on the post OP CT
% the result "elecModels" is a set of electrodes objects reconstructed from the CT, i.e. elecModels{1} is the first electrode found, elecModels{2} the second
elecModels = PaCER(niiCT); 

%% Transforming electrode objects from post OP CT space to T1 space
elecModels_transformed{1} = elecModels{1}.applyFSLTransform(FlirtMat, POSTOPCT_FILENAME, T1_FILENAME);
elecModels_transformed{2} = elecModels{2}.applyFSLTransform(FlirtMat, POSTOPCT_FILENAME, T1_FILENAME);

%% Plot MPR
figure('Name', 'Example 5.1');
%mpr = createSimpleMPRWorldCoordinates(niiT1); %scroll MPR planes by clicking a plane and turing the mouse wheel
%colormap gray;
hold on;

%% Plot electrodes
elecModels_transformed{1}.plot3D;
elecModels_transformed{2}.plot3D;

%% Plot Atlas Structures
atlas.plot3D;
atlas.setLabelColorByName(atlas.presentLabelNames, rgb('Goldenrod'));
atlas.setLabelColorByName({'STN_L', 'STN_R'}, rgb('LimeGreen'));

%% **** Continuing from EXAMPLE_5 here: ***

%% loading mrtrix tractography data (tensor-based tracking of CST fibers)
plotFraction = 0.1; % plot only XX share the fibres for performance reasons

mrtrix_CST_lh = read_mrtrix_tracks('native_nifti_data/TRACKS_CST_lh.t1space.tck');
mrtrix_CST_rh = read_mrtrix_tracks('native_nifti_data/TRACKS_CST_rh.t1space.tck');

fibSet1 = FiberSet(uniformRndSample(mrtrix_CST_lh.data, plotFraction));
fibSet2 = FiberSet(uniformRndSample(mrtrix_CST_rh.data, plotFraction));

fibSet1.plot3D;
fibSet2.plot3D;
view(30,10);

%% Add a Mädler/Coenen VTA projection from elec2 to fibset1
elecModels_transformed{2}.activeContact = 2;
vta = SimpleFieldSim3D(elecModels_transformed{2});
vta.impedance = 900;
vta.voltage = 9.6;
fibSet1.addSimpleFieldSim3D(vta);

%% Get a Config Windows
vta.getConfigWindow();
