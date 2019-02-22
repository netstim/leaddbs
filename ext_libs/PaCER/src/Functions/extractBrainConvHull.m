%% Extracts the convex hull of the brain out of a CT image
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dept. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2015 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu
function [convHullBrainMask, roughBrainMask] = extractBrainConvHull(niiCT)
%% force caching
niiCT.isToBeCached = 1;
if(~niiCT.isLoaded)
    niiCT.load();
end

%% Constants
LOWER_CT_BRAIN_THRESHOLD = 15; % [Hounsfield]
UPPER_CT_BRAIN_THRESHOLD = 60; % [Hounsfield], was: 60

% Check if the CT is in "standard" range [-1024 4096], if not 
% assume a 1024 offset was added
if(min(niiCT.img(:)) >= 0)
    LOWER_CT_BRAIN_THRESHOLD = LOWER_CT_BRAIN_THRESHOLD + 1024;
    UPPER_CT_BRAIN_THRESHOLD = UPPER_CT_BRAIN_THRESHOLD + 1024;
end
    
%% Algorithm
% downsample for performance improvment if needed
%disp('Downsampling is OFF. If this is to slow for you enable downsampling to get quicker (and dirtier) brain mask results. ');
ctIso = niiCT.img; % = downsampleImage(niiCT.img, niiCT.voxsize, niiCT.voxdim);

% threshold
ctIsoMedFilt = medfilt3(ctIso); 
% TODO in very rare cases (e.g. with neck coverage) the media filter might induce a connected componed 
% containing brain and some skin, which is super bad, however in most cases
% it is super good and makes most of the following processing unnessesary
% ;-)
% in principle this should be curable by using a distance transform and cutting the part that are weakly connect via
% such a weird path (having very long distnaces algong positive voxeles to
% the CoG) i.e. we end up with a fast marching approach

threImg = (ctIsoMedFilt > LOWER_CT_BRAIN_THRESHOLD & ctIsoMedFilt < UPPER_CT_BRAIN_THRESHOLD);
%threImg = reshape(threImg,niiCT.voxdim');
 [xx,yy,zz] = ndgrid(-2:2);
structEle = strel('sphere', ceil(2 / max(niiCT.voxsize)));  %sqrt(xx.^2 + yy.^2 + zz.^2) <= 2.5 / sqrt(max(niiCT.voxsize));
structEleSmall = sqrt(xx.^2 + yy.^2 + zz.^2) <= 1; %  CHECK THIS! (000)
%structEle = strel('disk', 2);
%G = gpuArray(structure); % if we have CUDA..

% morphology
%threImg = medfilt3(threImg);
morphImg = imopen(threImg,structEle);
morphFraction = sum(morphImg(:)) / numel(morphImg);
disp(['MorphFraction:' num2str(morphFraction)]);
if(morphFraction < 0.15 || morphFraction > 0.3)
    warning('Uncommon fraction of CT data in threshold range (15-60 HU). Trying to compensate. Make sure to use "soft tissue" reconstruction filters for the CT (e.g. J30 kernel) if this fails. ')
    %threImg = medfilt3(threImg);
    morphImg = imclose(threImg,structEle); %maybe we have a super noise brain tissue image, thus to CLOSE instead of open first
    morphImg = imerode(morphImg,structEle); % due to the closing all masks will enlarge, i.e. the brain mask might cover the ckull, thus we have to shrink it at the very end again
end 

clear threImg;

% largest connected component is brain % DO 2.5 instead of 3D for increased
% robustness ???
cc = bwconncomp(morphImg);
ccProps = regionprops(cc, 'Area', 'PixelIdxList');
[~, idx] = sort([ccProps.Area]);
brainMaskIdxs = ccProps(idx(end)).PixelIdxList;
roughBrainMask = false(size(morphImg));
%roughBrainMask(brainMaskIdxs) = true; % 3 D
roughBrainMask = largestConnCompSliceWise(morphImg); %2.5D
clear morphImg;

% axial slice wise convex hull (as the support of binary 3d conv hulls in
% matlab is bad and we don't need 100% accuracy for our task)
% and filling of rough mask
convHullBrainMask = false(size(roughBrainMask));
maskedCT = ctIsoMedFilt;
maskedCT(~roughBrainMask) = NaN; %threshold the roughBrainMask again to make sure no skull is contained due to morphology
roughBrainMask = (maskedCT > LOWER_CT_BRAIN_THRESHOLD & maskedCT < UPPER_CT_BRAIN_THRESHOLD); %TODO check if removable
for i=1:size(roughBrainMask,3)
    convHullBrainMask(:,:,i) = bwconvhull(roughBrainMask(:,:,i));
    roughBrainMask(:,:,i) = imfill(roughBrainMask(:,:,i), 'holes');
end

convHullBrainMask = imerode(convHullBrainMask,structEle);
convHullBrainMask = imerode(convHullBrainMask,structEle);

roughBrainMask = imerode(roughBrainMask,structEleSmall);

% convHullBrainMask = logical(upsampleImage(convHullBrainMask, niiCT.voxsize, niiCT.voxdim));
% roughBrainMask = logical(upsampleImage(roughBrainMask, niiCT.voxsize, niiCT.voxdim));

% 3D conv hull
% % [i,j,k]  = ind2sub(size(brainMask), find(brainMask));
% % 
% % i = i + ccProps(idx(end)).BoundingBox(1) * niiCT.voxsize(1);
% % j = j + ccProps(idx(end)).BoundingBox(2) * niiCT.voxsize(2);
% % k = k + ccProps(idx(end)).BoundingBox(3) * niiCT.voxsize(3);
% % 
% % [K,V]=  convhull(i,j,k);
% % trisurf(K, i, j, k)


end