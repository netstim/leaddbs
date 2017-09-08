%% electrodePointCloudModelEstimate - extract skeleton points from given pointcloud and subsequently fit a polynomial
%
% PARAMETERS
% varargin{1} = pixelValues list got from reference image of the  segmentation. When given used for intensity based
%               weightig of cenroid  calculation.
% varargin{2} = reverse Z direction (for special cases, e.g. certain phantom studies)
%
% Andreas Husch
% Centre Hospitalier de Luxembourg / Luxembourg Centre for Systems
% Biomedicine, University of Luxembourg
% 2014  - 2017
% mail@andreashusch.de

function [r3polynomial, tPerMm, skeleton, totalLengthMm] = electrodePointCloudModelEstimate(elecPointCloudMm, varargin)
revDir = false;
INTERNAL_DEGREE = 8; % fixed, determined as sufficient by AIC analysis

if(nargin < 2)
    USE_REF_IMAGE_WEIGHTING=false;
    warning('No Reference image for intensity weighting given! Accuracy is thus limited to voxel size!');
    pixelValues = [];
else
    USE_REF_IMAGE_WEIGHTING=true;
    pixelValues = varargin{1};
    if(nargin == 3)
        revDir = varargin{2};
    end
end

%% Axial centers based skeletonization
zPlanes = unique(elecPointCloudMm(:,3));
tol = 0;
if ~(length(zPlanes) < length(elecPointCloudMm))
    warning('CT planes in Z direction are not exactly aligned. Trying with 0.1 mm tolerance')
    tol = 0.1;
    zPlanes = uniquetol(elecPointCloudMm(:,3), tol / max(abs(elecPointCloudMm(:,3)))); % trying with tolerance
    
end

assert(length(zPlanes) < length(elecPointCloudMm), 'Couln''t find CT planes in z direction. Check that the CT scan was not acquired oblique!');

skeleton = [];%skeleton = NaN(length(zPlanes),3); length is unknown a-priori because of possible "non planes"
sumInPlane = []; %sumInPlane = NaN(length(zPlanes),1);

for i=1:length(zPlanes)
    inPlanePoints = elecPointCloudMm( abs(elecPointCloudMm(:,3) - zPlanes(i)) <= tol ,:);
    if(size(inPlanePoints,1) > 1)
        if(USE_REF_IMAGE_WEIGHTING)
            inPlaneIntensities = single(pixelValues(abs(elecPointCloudMm(:,3) - zPlanes(i)) <= tol)); % pixelValues MUST be same order than elecPointCloudMm!     

            skeleton(end+1,:) = inPlanePoints' * inPlaneIntensities / sum(inPlaneIntensities)'; %#ok<AGROW>
            sumInPlane(end+1) = sum(inPlaneIntensities); %#ok<AGROW> 
        else
            skeleton(end+1,:) = mean(inPlanePoints); %#ok<AGROW>
        end
    %else
        % ignore pseudo slices with just one plane, so do nothing 
    end
end

%% Filter Skeleton for valid Points
% see bar(sumInPlane) to get a feeling
filter = sumInPlane < (median(sumInPlane) / 1.5);
if(sum(filter) > 0)
    disp('Applied axial skeleton filter because of low intensity planes')
    skeleton = skeleton(~filter,:);
end

if(isequal(skeleton(1,:) , [0 0 0]))
    error('Empty skeleton. Was the CT image aquired in axial flow?')
end

%% Approximate parameterized polynomial  ([x y z] = f(t))
if(length(skeleton) < INTERNAL_DEGREE + 1)
    warning(['electrodePointCloudModelEstimate: less data points ' num2str(length(skeleton)) ...
        'than internal poly. degree (' num2str(INTERNAL_DEGREE) '). Lowering degree but take care']);
    INTERNAL_DEGREE = length(skeleton) - 1;
end
if(revDir)
    [r3polynomial, tPerMm] = fitParamPolyToSkeleton(flipud(skeleton), INTERNAL_DEGREE); % note degree (8)
else
    [r3polynomial, tPerMm] = fitParamPolyToSkeleton(skeleton, INTERNAL_DEGREE);  % note degree (8)  
end
totalLengthMm = polyArcLength3(r3polynomial, 0, 1); % for backwards compatibility
end


