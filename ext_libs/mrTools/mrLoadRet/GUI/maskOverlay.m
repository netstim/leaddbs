%
% maskOverlay.m
%     $Id$
%   usage: [maskOverlay, overlayData] = maskOverlay(thisView,<overlayList>,<scanList>,<boxInfo>)
%      by: jb, taken out of computeOverlay.m
%    date: 11/05/2010
% purpose: constructing a logical mask excluding clipped values based on clip values of each/all overlays in an analysis/scan
%   input:  - overlayList (default: all overlays in analysis) specifies which masks and overlays (and their associated alpha masks and overlays) must be output
%           - scanList (default: all scans in group): overlays are retrieved from the view for these scans
%           - boxInfo: this is a structure that specifies the transformation from a scan space to a slice/3D subset in another space
%               the outputs are then in this space
%
%           If clipAcrossOverlays is on, the mask is identical for all overlays (in a scan) and
%             is the intersection of all masks of all overlays (repeated so that maskOverlay is the same size as overlayData)
%             otherwise, one mask is output per overlay (on the fourth dimension) of overlayList

function [maskOverlays, overlays, overlayCoords] = maskOverlay(thisView,overlayList,scanList,boxInfo)

clipAcrossOverlays=viewGet(thisView,'clipAcrossOverlays');
nOverlaysInAnalysis = viewGet(thisView,'numberofOverlays');
analysisNum = viewGet(thisView,'currentAnalysis');

if ieNotDefined('scanList')
  scanList = 1:viewGet(thisView,'nScans');
end
if ieNotDefined('overlayList')
  overlayList=1:nOverlaysInAnalysis;
end

maskOverlays = cell(1,length(scanList));
overlays = cell(1,length(scanList));


if clipAcrossOverlays
  overlaysToGet = 1:nOverlaysInAnalysis;
else
  overlaysToGet = overlayList;
  overlayList = 1:length(overlaysToGet);
end

if ~ieNotDefined('boxInfo')
  boxInfo.interpMethod = mrGetPref('interpMethod');
  if isempty(boxInfo.interpMethod)
    boxInfo.interpMethod = 'linear';
  end
  boxInfo.interpExtrapVal = NaN;
  overlayCoords = cell(1,length(scanList));
end

cScan = 0;
overlayData = cell(1,length(scanList));
for iScan = scanList
  cScan = cScan+1;
  if ~ieNotDefined('boxInfo')
    % pull out overlays for this slice
    [overlayData{cScan},overlayCoords{cScan}] = ...
        getOverlayBox(thisView,iScan,boxInfo,analysisNum,overlaysToGet);
%     %put slices on 4th dimensions
%     overlayData{cScan} = permute(overlayData{cScan},[1 2 4 3]);
  else
    scanDims = viewGet(thisView,'scandims');
    overlayData{cScan}=NaN([scanDims nOverlaysInAnalysis]);
    cOverlay=0;
    for iOverlay = overlaysToGet
      cOverlay = cOverlay+1;
      thisOverlayData = viewGet(thisView,'overlayData',iScan,iOverlay);
      if ~isempty(thisOverlayData)
        overlayData{cScan}(:,:,:,cOverlay) = thisOverlayData;
      else
        overlayData{cScan}(:,:,:,cOverlay) = NaN;
      end
    end
  end
end

for iScan = 1:length(scanList)
  %we assume overlay-specific masks are all true for empty overlays 
  %(this is so that they (partially) empty overlays not mask everything when clipAcrossOverlays is on)
  maskOverlayData = true(size(overlayData{iScan}));
  
  cOverlay=0;
  for iOverlay = overlaysToGet
    cOverlay = cOverlay+1;
    thisOverlayData = overlayData{iScan}(:,:,:,cOverlay);
    
    if ~isempty(thisOverlayData) && ~all(isnan(thisOverlayData(:)))  %if there is some data in the overlay
      clip = viewGet(thisView,'overlayClip',iOverlay);
      
      if diff(clip) > 0 % Find defined pixels that are within clip
        maskOverlayData(:,:,:,cOverlay) = ((thisOverlayData >= clip(1) & thisOverlayData <= clip(2))) | isnan(thisOverlayData);
      elseif diff(clip) < 0 % Find defined pixels that are outside clip
        maskOverlayData(:,:,:,cOverlay) = (thisOverlayData >= clip(1) | thisOverlayData <= clip(2)) | isnan(thisOverlayData) ;
      else
        maskOverlayData(:,:,:,cOverlay) = false(size(thisOverlayData)) | isnan(thisOverlayData);
      end
    end
  end
  if clipAcrossOverlays % intersect all the overlay if clipAcross Overlays is on
    maskOverlayData=repmat(all(maskOverlayData,4),[1 1 1 size(maskOverlayData,4)]);
  end
  maskOverlays{iScan} = true([size(thisOverlayData,1) size(thisOverlayData,2) size(thisOverlayData,3) length(overlayList)]);
  maskOverlays{iScan}(:,:,:,logical(overlayList)) = maskOverlayData(:,:,:,overlayList(logical(overlayList)));
  overlays{iScan} = NaN([size(thisOverlayData,1) size(thisOverlayData,2) size(thisOverlayData,3) length(overlayList)],mrGetPref('defaultPrecision'));
  overlays{iScan}(:,:,:,logical(overlayList)) = overlayData{iScan}(:,:,:,overlayList(logical(overlayList)));
  
  %all non-defined voxels in the overlay must now become false in the mask
  maskOverlays{iScan}(isnan(overlays{iScan}))=false;
end



%%%%%%%%%%%%%%%%%%%%%%%
%   getOverlayBox   %
%%%%%%%%%%%%%%%%%%%%%%%
function [overlayImages,overlayCoords] = ...
  getOverlayBox(thisView,scanNum,boxInfo,analysisNum,overlayList)
%
% getOverlayBox: extracts overlays 3D subset and corresponding coordinates

overlayCoords = [];
overlayCoordsHomogeneous = [];
overlayImages = [];

base2overlay = boxInfo.base2overlay;
baseCoordsHomogeneous = boxInfo.baseCoordsHomogeneous;
baseDims = boxInfo.baseDims;
interpMethod = boxInfo.interpMethod;
interpExtrapVal = boxInfo.interpExtrapVal;

%find unique overlays in the list
[uniqueOverlays,dump,uniqueOverlaysIndices] = unique(overlayList);
nOverlays = length(uniqueOverlays);


interpFnctn = viewGet(thisView,'overlayInterpFunction',analysisNum);
% Transform base coords corresponding to this slice/image to overlays
% coordinate frame.
if ~isempty(base2overlay) & ~isempty(baseCoordsHomogeneous) 
  % Transform coordinates
  if size(baseCoordsHomogeneous,3)>1%if it is a flat map with more than one depth
    corticalDepth = viewGet(thisView,'corticalDepth');
    corticalDepthBins = viewGet(thisView,'corticalDepthBins');
    corticalDepths = 0:1/(corticalDepthBins-1):1;
    slices = corticalDepths>=corticalDepth(1)-eps & corticalDepths<=corticalDepth(end)+eps; %here I added eps to account for round-off erros
    baseCoordsHomogeneous = baseCoordsHomogeneous(:,:,slices);
    nDepths = nnz(slices);
  else
    nDepths=1;
  end
  overlayCoords = (base2overlay * reshape(baseCoordsHomogeneous,4,prod(baseDims)*nDepths))';
  overlayCoords = overlayCoords(:,1:3);
end


% Extract overlayImages

%if interpolation method is 'nearest', rounding the base coordinates 
%before interpolation should not change the result
%then, using unique to exclude duplicate coordinates,
%interpolation can be performed on far less voxels (at least for surfaces, 
%or if the base has a higher resolution than the scan)
%and this can substantially reduce the time needed for interpolation
if ~isempty(overlayCoords) && strcmp(interpMethod,'nearest') 
  overlayCoords = round(overlayCoords);   
  [overlayCoords,dump,coordsIndex] = unique(overlayCoords,'rows'); 
end
  
if ~isempty(interpFnctn)
  overlayImages = feval(interpFnctn,thisView,scanNum,baseDims,...
    analysisNum,overlayCoords,interpMethod,interpExtrapVal,uniqueOverlays);
else
  overlayImages = zeros([size(overlayCoords,1) nOverlays]);
  cOverlay=0;
  for iOverlay = uniqueOverlays
    cOverlay = cOverlay+1;
    if iOverlay
      overlayData = viewGet(thisView,'overlayData',scanNum,iOverlay,analysisNum);
      if ~isempty(overlayData) & ~isempty(overlayCoords)
        % Extract the slice
        overlayImages(:,cOverlay) = interp3(overlayData,...
          overlayCoords(:,2),overlayCoords(:,1),overlayCoords(:,3),...
          interpMethod,interpExtrapVal);
      else
        overlayImages(:,cOverlay) = NaN;
      end
    else
      overlayImages(:,cOverlay) = NaN;
    end
  end
end

if ~isempty(overlayImages)
  if strcmp(interpMethod,'nearest')
    overlayImages = reshape(overlayImages(coordsIndex,:,:),[baseDims(1) baseDims(2) baseDims(3)*nDepths nOverlays]);
    overlayCoords = reshape(overlayCoords(coordsIndex,:,:),[baseDims(1) baseDims(2) baseDims(3)*nDepths 3]);
  else
    overlayImages = reshape(overlayImages,[baseDims(1) baseDims(2) baseDims(3)*nDepths nOverlays]);
    overlayCoords = reshape(overlayCoords,[baseDims(1) baseDims(2) baseDims(3)*nDepths 3 ]);
  end
  if nDepths>1
    %we take the average depth of all displayed cortical depths as the actual cortical depth location
    overlayCoords = mean(overlayCoords,3);
  end
  overlayImages = overlayImages(:,:,:,uniqueOverlaysIndices);
end


%%%%%%%%%%%%%%%%%%%%%%%
%    corAnalInterp   %%
%%%%%%%%%%%%%%%%%%%%%%%
function overlayImages = corAnalInterp(thisView,scanNum,...
  baseDims,analysisNum,overlayCoords,interpMethod,interpExtrapVal,overlayList)
%
% corAnalInterp: special case for corAnal. Need to treat amp/ph as complex
% valued.

if ieNotDefined('overlayList')
  overlayList = 1:viewGet(thisView,'numberofoverlays',analysisNum);
end

% Initialize
overlayImages = zeros(size(overlayCoords,1), length(overlayList));

cOverlay=0;
for iOverlay = overlayList
  cOverlay = cOverlay+1;
  if iOverlay
    overlayData = viewGet(thisView,'overlayData',scanNum,iOverlay,analysisNum);
    %check if overlay is an amplitude or a phase
    % if it s the case, load both and interpolate the complex value
    % we assume that phase overlays always follow amplitude overlays
    % Doing a bit more of book keeping, we could avoid interpolating twice 
    % the phase and amplitude overlays, 
    % but I don't think it's worth the trouble
    overlayType = viewGet(thisView,'overlayType',iOverlay,analysisNum);
    switch(overlayType)
      case 'amp'
        if ~strcmp(viewGet(thisView,'overlayType',iOverlay+1,analysisNum),'ph')
          mrWarnDlg(['(corAnalInterp) Could not find phase associated to overlay' viewGet(thisView,'overlayName',iOverlay,analysisNum) ', using non-complex interpolation']);
        else
          amplitude = overlayData;
          phase = viewGet(thisView,'overlayData',scanNum,iOverlay+1,analysisNum);
          overlayData = amplitude .* exp(1i*phase);
        end
      case 'ph'
        phase = overlayData;
        if ~strcmp(viewGet(thisView,'overlayType',iOverlay-1,analysisNum),'amp')
          mrWarnDlg(['(corAnalInterp) Could not find amplitude associated to overlay ' viewGet(thisView,'overlayName',iOverlay,analysisNum) ', using non-weighted phase interpolation']);
          overlayData = exp(1i*phase);
        else
          amplitude = viewGet(thisView,'overlayData',scanNum,iOverlay-1,analysisNum);
          overlayData = amplitude .* exp(1i*phase);
        end
    end

    if ~isempty(overlayData) && ~isempty(overlayCoords)
      % Extract the slice
      overlayImages(:,cOverlay) = interp3(overlayData,...
        overlayCoords(:,2),overlayCoords(:,1),overlayCoords(:,3),...
        interpMethod,interpExtrapVal);
      switch(overlayType)
        case 'amp'
          overlayImages(:,cOverlay) = abs(overlayImages(:,cOverlay));

        case 'ph'
          ang = angle(overlayImages(:,cOverlay));
          ang(ang < 0) = ang(ang < 0)+2*pi;
          overlayImages(:,cOverlay) = ang;
      end
    else
      overlayImages(:,cOverlay) = NaN;
    end
    
  else
    overlayImages(:,cOverlay) = NaN;
  end
end

