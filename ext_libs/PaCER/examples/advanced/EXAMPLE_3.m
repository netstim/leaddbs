%% PaCER/EXAMPLE_3.m - Plan/Outcome Comparisions
%
% Andreas Husch
% Centre Hospitalier de Luxembourg (CHL) / Luxembourg Centre for Systems
% Biomedicine (LCSB), University of Luxembourg
% mail@andreashusch.de / husch.andreas@chl.lu
% 2017
%
%% Make sure all PaCER files (including subdirectories) are in your MATLAB path!
SETUP_PACER % SETUP_PACER has to be run *once* to setup the PaCER path, further runs are not required

%% Load post OP CT 
POSTOPCT_FILENAME = '8e7f749_CT_POSTOP_short.nii.gz'; %low res immediate post OP Scan (brain shift still present), co-registered to pre OP planning CT space
XML_PLAN = '8e7f749.xml';
niiCT = NiftiMod(POSTOPCT_FILENAME);

%% Run PACER
[elecModels, elecPointCloudsStruct] = PaCER(niiCT, 'medtronicXMLPlan', XML_PLAN, 'electrodeType', 'Medtronic 3389'); % by supplying the XML_PLAN pacer ...
% ...automatically associates reconstructed electrodes to planned
% trajectories, PaCER will issue a warning due to low SNR of the
% low-resolution post OP scan, however this is fine in this case as we are
% interested in the **trajectory** recontruction which is still accuracte

%% Plot MPR of CT
figure('name', 'Example 3');
mpr = createSimpleMPRWorldCoordinates(niiCT); %scroll MPR planes by clicking a plane and turing the mouse wheel
colormap gray;
camzoom(2);
hold on;

%% Plot electrodes
elecModels{1}.plot3D;
elecModels{2}.plot3D;

%% Turn sterotactic plan(s) into a plottable TestElectrodes object
% Note that the ben's gun **rotation** is approximate as it depends on the
% kinematics of the sterotactic frame / microdrive device
testElecs = {};
for i=1:2
    e1 = convertMedtronicCoordToLPI(elecPointCloudsStruct(i).associatedXmlDefinition.entry, niiCT);
    t1 = convertMedtronicCoordToLPI(elecPointCloudsStruct(i).associatedXmlDefinition.target, niiCT);
    
    e1world = niiCT.getNiftiWorldCoordinatesFromMatlabIdx(e1 ./ niiCT.voxsize);
    t1world = niiCT.getNiftiWorldCoordinatesFromMatlabIdx(t1 ./ niiCT.voxsize);
    
    testElecs{end+1} = TestElectrodes(e1world,t1world); %#ok<SAGROW>
    testElecs{end}.plot3D;
end
