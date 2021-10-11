%% PaCER - Precise and Convenient Electrode Reconstruction for Deep Brain Stimulation
%
% Andreas Husch
% Centre Hospitalier de Luxembourg (CHL), Dept. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine (LCSB)
% (c) 2016 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu

function [elecModels, elecPointCloudsStruct, intensityProfiles, skelSkelmms] = PaCER(niiCT, varargin)
%% Optional Arguments and Default Values
argParser = inputParser();
argParser.KeepUnmatched = true;

argParser.addOptional('finalDegree', 3, @(x)(isnumeric(x) && (x >= 1)));

argParser.addOptional('displayProfiles', false); % optional plot of intensity profiles
argParser.addOptional('displayMPR', false); % optional MPR plot of orthogonal oblique resampling along the trajceotry

argParser.addOptional('noMask', false); % for phantom studies where no brain is present in data
argParser.addParameter('brainMask', ''); % for manually providing a brain mask (binary segmentation image file path)

argParser.addOptional('reverseDir', false); % for special cases with I-S flip
argParser.addOptional('contactDetectionMethod', 'contactAreaCenter', @(x)(ismember(x, {'peak', 'peakWaveCenter', 'contactAreaCenter'}))); % default contactAreaCenter, if peak, automatic fallback to contactAreaCenter for "bad quality data"
argParser.addParameter('electrodeType', '', @(x)(ismember(x, {'', 'Medtronic 3387', 'Medtronic 3389', 'Boston Vercise Directional', 'Aleva Neurotherapeutics directSTIM Directed'}))); 

argParser.addOptional('medtronicXMLPlan', '', @(x)(ischar(x))); 
argParser.parse(varargin{:});
args = argParser.Results;

%% profe contactDetectionMethod contactAreaCenter if electrodeType is set manually
if(~isempty(args.electrodeType))
   args.contactDetectionMethod = 'contactAreaCenter'; 
end
%% Checks
assert(logical(license('test', 'image_toolbox')), 'It seems this system does not have the Image Processing Toolbox installed. PaCER requires the Image Processing Toolbox to continue.')

if(~isa(niiCT, 'NiftiMod') && ~isa(niiCT, 'NiftiModSPM') )
    disp('First parameter is not a nifti object. Intrepretating as filename and tring to load a nifti file with that name from disk...');
    niiCT = NiftiMod(niiCT);
end

if(max(niiCT.voxsize) > 1)
    warning('Slice thickness is greater than 1 mm! Independent contact detection is most likly not possible. Forcing contactAreaCenter based method.');
    args.contactDetectionMethod = 'contactAreaCenter';
elseif(max(niiCT.voxsize) > 0.7)
    warning('Slice thickness is greater than 0.7 mm! Independet contact detection might not work reliable in this case. However, for certain electrode types with large contacts spacings you might be lucky.');
end

%% Run Algorithm
disp(['===================== Processing ' niiCT.filepath ' =====================']);
disp(['Voxel size: ' num2str(niiCT.voxsize')]);

elecPointCloudsStruct = extractElectrodePointclouds(niiCT, varargin{:}); % preprocessing
elecModels = {};
intensityProfiles = {};
skelSkelmms = {};
for i=1:length(elecPointCloudsStruct) 
    disp('------------- Processing Electrode -------------');
    % Preprocessing and "1st pass" model
    initialR3polynomial = electrodePointCloudModelEstimate(elecPointCloudsStruct(i).pointCloudWorld , elecPointCloudsStruct(i).pixelValues, args.reverseDir); % internally always degree 8
    % Refitting ("2nd and 3rd pass")
    [elecModels{i}, intensityProfiles{i}, skelSkelmms{i}] = refitElec(initialR3polynomial,elecPointCloudsStruct(i).pointCloudWorld, elecPointCloudsStruct(i).pixelValues, args); %#ok<AGROW>
end
end