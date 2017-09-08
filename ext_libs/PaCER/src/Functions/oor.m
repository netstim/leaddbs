%% oor - optimal opliptimal oblique sampling along a given polynonimal
%
% Returns: improvedSkeleton - skeleton gained along the polynomial by
%                             resampling and intensity weithed centroid calculation
%          medIntensity     - median intensity profile along the polynoial
%
% Andreas Husch
% Centre Hospitalier de Luxembourg / Luxembourg Centre for Systems
% Biomedicine, University of Luxembourg
% 2016  - 2017
% mail@andreashusch.de
%
% That code needs some cleanup :-) But it works ;-)
function [improvedSkeleton, medIntensity,orthIntensVol, orthSamplePointsVol, skelScaleMm] = orr(r3polyToUse, STEP_SIZE, XGrid, YGrid, interpolationF)
SND_THRESH=1500;

arcLength = polyArcLength3(r3polyToUse, 0, 1);
oneMmEqivStep = 1 / arcLength;
lookahead = 3 * oneMmEqivStep; % 2 mm lookahead

regX = r3polyToUse(:,1);
regY = r3polyToUse(:,2);
regZ = r3polyToUse(:,3);
orthogonalSamplePoints = []; %#ok<NASGU>
improvedSkeleton = [];
avgIntensity = [];
medIntensity = [];
sumIntensity = [];
ind=1;
% cf. invPolyArcLength3(r3polyToUse, 0:Z_RESOLUTION:arcLength)) (expensive,
%  DOUBLE CHECK invPolyArcLength3 vs STEP_SIZE     CHECK RES: For slightly
%  beended elecs it seems accurate below 10 micron
% diff(..) is constant for low resolution / limited bending, i.e. the
% following is in this case a very good approximation
evalAtT =-lookahead:STEP_SIZE:1; % run a bit further than the (currently believed) tip (look-a-head)
for evalAt = evalAtT  
    x_d = polyval(polyder((regX)), evalAt);
    y_d = polyval(polyder((regY)), evalAt);
    z_d = polyval(polyder((regZ)), evalAt);
    
    currentPoint = polyval3(r3polyToUse,evalAt)';
    direction = [x_d y_d z_d];
    directionNorm = direction / norm(direction);
    
    ortho1 = cross(directionNorm', [0 1 0]');
    ortho1 = ortho1 / norm(ortho1);
    ortho2 = cross(directionNorm', ortho1);
    ortho2 = ortho2 / norm(ortho2);
    orthogonalSamplePoints = bsxfun(@plus, currentPoint, ortho1 * XGrid(:)' + ortho2 * YGrid(:)'); % XGrid, YGrid in mm
    orthSamplePointsVol(:,:,ind) = orthogonalSamplePoints; %#ok<AGROW>
    intensities = interpolationF(orthogonalSamplePoints');
    intensitiesNanZero = intensities;
    intensitiesNanZero(isnan(intensitiesNanZero)) = 0;
    intensitiesNanZero(intensities<SND_THRESH) = 0;
  
    % determine new skel point after rethresholding
    skelPoint = orthogonalSamplePoints * intensitiesNanZero / sum(intensitiesNanZero);
    if(any(isnan(skelPoint)))
        evalAtT(ind) = [];
        continue;
    end
    avgIntensity(ind) = nanmean(intensities); %#ok<AGROW>
    sumIntensity(ind) = nansum(intensitiesNanZero); %#ok<AGROW>
    medIntensity(ind) = nanmedian(intensities); %#ok<AGROW>
    
    improvedSkeleton(ind,:) = skelPoint; %#ok<AGROW>
    intensityMap = reshape(intensitiesNanZero,size(XGrid));

    orthIntensVol(:,:,ind) = intensityMap; %#ok<AGROW>
    
    ind=ind+1;    
end
lowerLimits = zeros(length(evalAtT),1);
upperLimits = evalAtT;
lowerLimits(upperLimits < 0) = upperLimits(upperLimits < 0);
upperLimits(upperLimits < 0) = 0;
skelScaleMm = polyArcLength3(r3polyToUse, lowerLimits, upperLimits); 
skelScaleMm ( lowerLimits < 0) = -skelScaleMm ( lowerLimits < 0) ;
end