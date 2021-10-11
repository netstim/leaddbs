%% refitElec - refits a previously extracted corse electrode model using
%               ptimal oblique samplingto increase accuracy and enable 
%               automatic contact localisation. Intital refitting
%               ("2nd pass" is followed by a final "3rd pass" to
%               calibrate the point of origin / zero according to the detected contacts)
%
% Params: initialPoly     - Dx3 double, R3 ploynomial coeffecients, 
%         pointCloudWorld - Nx3 double, point cloud of sourounding voxels in respective coordinate system
%         voxelValues     - Nx1 double, voxels intensitiy values (1:1 to  pointCloudWorld)
%         varargin        - options cell
%
% Returns: refitReZeroedElecMod - polynomialElectrodeModel Object with the
%                                  refitted and rezeored elec model,
%                                  including reference to elecInfo etc.
%
% Andreas Husch
% Centre Hospitalier de Luxembourg / Luxembourg Centre for Systems
% Biomedicine, University of Luxembourg
% 2016  - 2017
% mail@andreashusch.de

function [refitReZeroedElecMod, filteredIntensity, skelScaleMm] = refitElec(initialPoly, pointCloudWorld, voxelValues, varargin)
%% CONSTANTS FOR OOR / CONTACT DETECTION
XY_RESOLUTION = 0.1;
Z_RESOLUTION = 0.025;
LIMIT_CONTACT_SEARCH_MM = 20; % limit contact search to first 20 mm

%% Options and defaults
argParser = inputParser();
argParser.KeepUnmatched = true;

argParser.addParameter('finalDegree', 3);
argParser.addParameter('contactDetectionMethod', 'peak', @(x)(ismember(x, {'peak', 'peakWaveCenter', 'contactAreaCenter'}))); % optional MPR plot of orthogonal oblique resampling along the trajceotry

argParser.addParameter('displayProfiles', false); % optional plot of intensity profiles
argParser.addParameter('displayMPR', false); % optional MPR plot of orthogonal oblique resampling along the trajceotry

argParser.addParameter('electrodeType', '', @(x)(ismember(x, {'', 'Medtronic 3387', 'Medtronic 3389', 'Boston Vercise Directional', 'Aleva Neurotherapeutics directSTIM Directed'}))); 

argParser.parse(varargin{:});
args = argParser.Results;

FINAL_DEGREE = args.finalDegree;
DISPLAY_PROFILES = args.displayProfiles;
DISPLAY_MPR = args.displayMPR;

%% Algorithm
interpolationF = scatteredInterpolant(pointCloudWorld, double(voxelValues), 'linear', 'none'); %TODO check cubic!
totalLengthMm  = polyArcLength3(initialPoly, 0, 1);

[XGrid,YGrid] = meshgrid(-1.5:XY_RESOLUTION:1.5,-1.5:XY_RESOLUTION:1.5); % 4 mm  around initial skel with 0.1 mm resolution
oneMmEqivStep = 1 / totalLengthMm; % approx.
STEP_SIZE = Z_RESOLUTION * oneMmEqivStep; % mm * toMm

%% 2nd Pass
skeleton2nd = oor(initialPoly, STEP_SIZE, XGrid, YGrid, interpolationF);
refittedR3Poly2nd = fitParamPolyToSkeleton(skeleton2nd,8);  %<== "INTERNAL" DEGREE fittet including lookahead!
% this additional call is only getting the skelScaleMm and medIntensity for the 2nd pass refitted Poly.
[skeleton3rd, medIntensity, orthIntensVol, ~, skelScaleMm] = oor(refittedR3Poly2nd, STEP_SIZE, XGrid, YGrid, interpolationF);
%[~, ~, medIntensity, ~,~,orthIntensVol, ~, skelScaleMm] = getOrthogonalCogSkel(refittedR3Poly2nd, STEP_SIZE, XGrid, YGrid, interpolationF);

disp(['1st Pass Electrode Length within Brain Convex Hull: ' num2str(polyArcLength3(initialPoly, 0, 1)) 'mm']);
disp(['2nd Pass Electrode Length within Brain Convex Hull: ' num2str(polyArcLength3(refittedR3Poly2nd, 0, 1)) 'mm']);

%figure, imagesc(intensityMap);
if (DISPLAY_MPR)
    figure, MPR(gca, orthIntensVol), daspect(1./[XY_RESOLUTION XY_RESOLUTION Z_RESOLUTION]) % make sure the daspect matches the x,y,z resolution ratios
end
%save
%figure, plot(scaleMm, sumIntensity)
%figure, plot(scaleMm, avgIntensity) % quite different to sum which is unexpected!?



%% Find 1D Intensity Peaks (i.e. electrode contacts / X-Ray markers)
% not that polynomial orthogonal skeletonizsation might run in "negative depth" (look-a-head to
% make sure to find tip!) however, the interpolation is precisly only on
% domain 0 1!
filterWidth = (0.25 / Z_RESOLUTION) + 1; % [mm] => [samples]

filteredIntensity = filtfilt(ones(1,filterWidth),filterWidth,medIntensity); % Since Matlab 2016a: movmean(medIntensity, filterWidth); % check med vs avg :)
filterIdxs = find(skelScaleMm <= LIMIT_CONTACT_SEARCH_MM); % restrict to first 20mm to save unnessecary computations
[peakLocs, peakWaveCenters, peakValues, threshIntensityProfile, threshold, contactAreaCenter, contactAreaWidth, xrayMarkerAreaCenter, xrayMarkerAreaWidth] = getIntensityPeaks(filteredIntensity, skelScaleMm, filterIdxs);

%% Decide on contact detection method: TODO refactor
if(strcmp(args.contactDetectionMethod, 'peakWaveCenters'))
    contactPositions = peakWaveCenters;
else %strcmp(args.contactDetectionMethod, 'peaks'))
    contactPositions = peakLocs;
end

if(DISPLAY_PROFILES)
    plotIntensityProfileAndPeaks(filteredIntensity, skelScaleMm);
end

disp(['refitElec: selected contactDetectionMethod is ' args.contactDetectionMethod]);
try
    [electrodeInfo, dataModelPeakRMS] = determineElectrodeType(contactPositions); % TODO fallback if unknown!
catch
    disp('Falling back to contactDectionMethod = "contactAreaCenter"');
    args.contactDetectionMethod = 'contactAreaCenter';
end

%% Find Zero / Point of Origin TODO: refactor, DRY with previous cell!
useDetectedContactPositions = 1;
if(length(contactPositions) < 4 || strcmp(args.contactDetectionMethod, 'contactAreaCenter'))
    %explicit request for contactAreaCenter and manual electroe type or fallback to due signal quality
    useDetectedContactPositions = 0;
    if(length(contactPositions) < 4)
        warning('Could NOT detect independent electrode contacts. Check image quality. ');
    end
    electrodeGeometries = load('electrodeGeometries.mat');
    electrodeGeometries = electrodeGeometries.electrodeGeometries;
    
    if(isempty(args.electrodeType)) % no electrode type given TODO: DRY!
        disp('Trying to estimate electrode type by simple contactAreaCenter. ');
        warning('No electrode specification given! Set electrodeType option! Trying to estimate type by contactAreaWidth only which might be wrong!');
        if(contactAreaWidth < 10.5)
            disp('Assuming Medtronic 3389 or Boston Scientific Vercise Directional. Setting 3389.');
            electrodeInfo  = electrodeGeometries(1);
        else
            disp('Assuming Medtronic 3387. Setting 3387.');
            electrodeInfo  = electrodeGeometries(2);
        end   
    else
        disp(['Setting user specified electrode type (' args.electrodeType ')']);
        [flag, idx] = ismember(args.electrodeType, {electrodeGeometries.string});
        electrodeInfo  = electrodeGeometries(idx);
        if(~flag)
            error('Unknown electrode type given');
        end
    end
    zeroT = invPolyArcLength3(refittedR3Poly2nd, contactAreaCenter-mean(electrodeInfo.ringContactCentersMm)); % calibrate zero

%     refittedR3PolyReZeroed = fitParamPolyToSkeleton(polyval3(refittedR3Poly2nd, linspace(zeroT,1,totalLengthMm  / XY_RESOLUTION)'),FINAL_DEGREE); % <=== LOWER DEGREE IS BETTER (FOR TIP!!) (TRADE OF)
%     
%     refitReZeroedElecMod = PolynomialElectrodeModel(refittedR3PolyReZeroed);
%     refitReZeroedElecMod.useDetectedContactPositions = 0; %
else % Indivdual Contact Based (i.e. detect electrode type automaticallay and use first contact)
    if(dataModelPeakRMS > 0.3) % TODO the MAX deviation might be a better measure than the RMS?
        disp('Switching to model based contact positions because of high RMS (Setting useDetectedContactPositions = 0).');
        useDetectedContactPositions = 0;
    end
    disp([electrodeInfo.string ' Electrode']);
    zeroT = invPolyArcLength3(refittedR3Poly2nd, contactPositions(1)-electrodeInfo.zeroToFirstPeakMm);
    refittedContactDistances = contactPositions - (contactPositions(1)-electrodeInfo.zeroToFirstPeakMm);

    %% Plot Profiles and Detected Peaks/Contacts
    %plot(skelScaleMm, medIntensity2ndThresh)
    %hold on,
    % plot(skelScaleMm, medIntensity)
    % hold on, plot(skelScaleMm, avgIntensity)
    %
    % %hold on, plot(skelScaleMm, medfilt1(avgIntensity,filterWidth))
    % hold on, plot(skelScaleMm, movmean(medIntensity, filterWidth))
    % hold on, plot(skelScaleMm, filtfilt(ones(1,filterWidth)./filterWidth,1,medIntensity))
    %
    % legend('Median Intens', 'Avg. Intens.', 'Filt. Med. Intens.', 'Zero-Phase Filt') % 'Median Intens. Re-Thresh.',
    % title('Electrode Trajectory Intensity Profiles')
    % ax.GridLineStyle = ':';
    % axis tight
    % grid on
    % figure('Name', 'Intensity Profile');
    % ax = gca;
    % xlabel('Length [mm]');
    % ylabel('Intensity Measure [HU]');
    % ylim([0 3200]);
    % xlim([0 20]) % show only first 20 mm
    % hold on;
    % hProfile = plot(skelScaleMm, filteredIntensity);
    % plot(peakLocs, peakValues + 0.02 * peakValues, 'v', 'MarkerFaceColor', hProfile.Color, 'MarkerEdgeColor', hProfile.Color);
    % ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
    % plot(skelScaleMm(filterIdxs), threshIntensityProfile(filterIdxs));
    % scatter(peakWaveCenters, repmat(threshold,1,4), 'filled');
    % grid on;
end

if(exist('detectElectrodeRotationBostonVercise.m', 'file')) % experemential lead rotation detection module is present
    if(strcmp(electrodeInfo.string, 'Boston Vercise Directional') || strcmp(args.electrodeType, 'Boston Vercise Directional'))
        try
            [rotationVector, markerPoint] = detectElectrodeRotationBostonVercise(skeleton2nd, skelScaleMm, xrayMarkerAreaCenter, xrayMarkerAreaWidth);
            electrodeInfo.rotationVector = rotationVector;
            electrodeInfo.rotationMarkerPoint = markerPoint;
        catch
            warning('Rotation detection failed');
        end
    end
end
%%% Refit Poly To Correct Zero Point (lower end of lowest electrode contact)
% Determine new "0 Point" ==def 1: lower end of the lowest electrode contact
% i.e. center of lowest electrode contact -0.75mm in proximal direction for
% 3389 / 3387, proximal direction is approximate!
if(FINAL_DEGREE == 1)
    elecEndT = invPolyArcLength3(refittedR3Poly2nd, LIMIT_CONTACT_SEARCH_MM); % fit only contact region but use original length in case of deg. 1
    refittedR3PolyTmp = fitParamPolyToSkeleton(polyval3(refittedR3Poly2nd, linspace(zeroT,elecEndT,totalLengthMm  / XY_RESOLUTION)'),FINAL_DEGREE); % <=== LOWER DEGREE IS BETTER (FOR TIP!!) (TRADE OF)
    refittedR3PolyReZeroed = fitParamPolyToSkeleton(polyval3(refittedR3PolyTmp, linspace(0,invPolyArcLength3(refittedR3PolyTmp,totalLengthMm),totalLengthMm  / XY_RESOLUTION)'),FINAL_DEGREE); % <=== LOWER DEGREE IS BETTER (FOR TIP!!) (TRADE OF)
else % normal case
    refittedR3PolyReZeroed = fitParamPolyToSkeleton(polyval3(refittedR3Poly2nd, linspace(zeroT,1,totalLengthMm  / XY_RESOLUTION)'),FINAL_DEGREE); % <=== LOWER DEGREE IS BETTER (FOR TIP!!) (TRADE OF)
end

disp(['Electrode Length within Brain Convex Hull after contact detection and Zero-Point calibration: ' num2str(polyArcLength3(refittedR3PolyReZeroed, 0, 1)) 'mm']);

refitReZeroedElecMod = PolynomialElectrodeModel(refittedR3PolyReZeroed, electrodeInfo);

if(strcmp(refitReZeroedElecMod.electrodeInfo.string,'Unkown Electrode Type') ||  useDetectedContactPositions == 0)
    refitReZeroedElecMod.useDetectedContactPositions = 0;
else
    refitReZeroedElecMod.useDetectedContactPositions = 1; 
    refitReZeroedElecMod.detectedContactPositions = refittedContactDistances(1:electrodeInfo.noRingContacts,:)';
end

