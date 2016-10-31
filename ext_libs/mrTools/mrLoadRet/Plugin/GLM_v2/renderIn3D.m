% renderIn3D - displays 3D rendering of currently loaded ROIs and Overlays
%
%        $Id$
%      usage: [  ] = renderRois3D(handle)
%         by: julien besle
%       date: 2011-08-09
%
%    purpose: displays 3D rendering of currently loaded ROIs and Overlays
%             (callback)


function renderIn3D(hObject, dump)
tic
handles = guidata(hObject);
viewNum = handles.viewNum;
thisView = viewGet(viewNum,'view');

scanNum = viewGet(thisView,'currentScan');
overlayList = viewGet(thisView,'currentOverlay');
roiList = viewGet(thisView,'visibleRois');
nRois = length(roiList);

if ~nRois
   mrWarnDlg('(renderRois3D) Please make at least one ROI visible');
   return;
end

% get the analysis structure
analysis = viewGet(thisView,'analysis');
currentAnalysis = viewGet(thisView,'currentAnalysis');
overlayAlpha = viewGet(thisView,'alpha');

%------------------------------------ Get some info on the base
baseNum = viewGet(thisView,'currentBase');
d.baseType = viewGet(thisView,'baseType',baseNum);
switch(d.baseType)
  case 0 %if the base is a volume
    baseDims = viewGet(thisView,'basedims',baseNum);
    overlayDims = viewGet(thisView,'scanDims');
    d.base2scan = viewGet(thisView,'base2scan');
    baseVoxelSize = viewGet(thisView,'baseVoxelSize');
  case 1 %if the base is a flat map
    baseCoordMap = viewGet(thisView,'baseCoordMap',baseNum);
    baseDims = viewGet(thisView,'basedims',baseNum);
    corticalDepthBins = viewGet(thisView,'corticalDepthBins');
    corticalDepths = viewGet(thisView,'corticaldepth');
    corticalDepthIndices = find((0:1/(corticalDepthBins-1):1)>=corticalDepths(1)-eps & (0:1/(corticalDepthBins-1):1)<=corticalDepths(end)+eps);
    baseDims(3) = length(corticalDepthIndices);
    rotateAngle = viewGet(thisView,'rotate');
    [baseImage,~,baseCoordsHomogeneous] = ...
          getBaseSlice(thisView,viewGet(thisView,'curslice'),viewGet(thisView,'baseSliceIndex',baseNum),rotateAngle,baseNum,d.baseType);
    baseCoordsHomogeneous = baseCoordsHomogeneous(:,:,corticalDepthIndices);
    baseCoordsHomogeneous = reshape(baseCoordsHomogeneous,4,[]);
    baseCoordsHomogeneous = round(baseCoordsHomogeneous);
    baseCoordsLinear = mrSub2ind(baseCoordMap.dims,baseCoordsHomogeneous(1,:),baseCoordsHomogeneous(2,:),baseCoordsHomogeneous(3,:));  %linearize the index
    d.base2scan = eye(4);
    overlayDims = baseDims;
    %compute the mean distance between inner and outer surface coordinates in mm
    meanThickness = mean(mean(sum((baseCoordMap.innerCoords - baseCoordMap.outerCoords).^2,4)));
    %in order to get an approximation of the scaling between axes
    flatResolution = 3; %we assume that flat maps have a resolution of 3 times the base voxel dimensions
    baseVoxelSize = [1/flatResolution 1/flatResolution meanThickness/corticalDepthBins];
  case 2
    d.base2scan = viewGet(thisView,'base2scan');
    overlayDims = viewGet(thisView,'scanDims');
    baseVoxelSize = viewGet(thisView,'baseVoxelSize');
end


%------------------------------------ Construct ROI objects ------------------------------------------------%
fprintf(1,'Constructing 3D ROIs...');
edgeX = [-.5 -.5; -.5 -.5; -.5 -.5; -.5 -.5; -.5 .5; -.5 .5; -.5 .5; -.5 .5; .5 .5; .5 .5; .5 .5; .5 .5];
edgeY = [-.5 .5; -.5 .5; -.5 -.5; .5 .5; -.5 -.5; -.5 -.5; .5 .5; .5 .5; -.5 .5; -.5 .5; -.5 -.5; .5 .5];
edgeZ = [-.5 -.5; .5 .5; -.5 .5; -.5 .5; -.5 -.5; .5 .5; .5 .5; -.5 -.5; -.5 -.5; .5 .5; -.5 .5; -.5 .5];
cubeEdgesCoords = [edgeX edgeY edgeZ];

allRoisCoords = [];
for iRoi = 1:nRois

  roi = viewGet(thisView, 'roi',roiList(iRoi));
  if d.baseType==1
    thisRoiCoords = getROICoordinates(thisView,roiList(iRoi),0,[],'baseNum',baseNum); %get the coordinates of the ROi in the canonical base space
    roiBaseCoordsLinear = mrSub2ind(baseCoordMap.dims,thisRoiCoords(1,:),thisRoiCoords(2,:),thisRoiCoords(3,:)); %linearize the index
    baseIndices = find(ismember(baseCoordsLinear,roiBaseCoordsLinear));
    [roiCoords{iRoi}(:,1),roiCoords{iRoi}(:,2),roiCoords{iRoi}(:,3)] = ind2sub(baseDims,baseIndices);
    roiCoords{iRoi} = roiCoords{iRoi}';
  else
    roiCoords{iRoi} = getROICoordinates(thisView,roiList(iRoi));
  end
  roiColor(iRoi,:) = color2RGB(roi.color);
  roiName{iRoi} = roi.name;
  
  if ~isempty(roiCoords{iRoi}) && d.baseType==0
        % we need to remove any coordinate that might fall outside the base anatomy (or do we ?)
        outside_voxels = roiCoords{iRoi}(1,:)<1 | roiCoords{iRoi}(1,:)>baseDims(1) |...
                             roiCoords{iRoi}(2,:)<1 | roiCoords{iRoi}(2,:)>baseDims(2) |...
                             roiCoords{iRoi}(3,:)<1 | roiCoords{iRoi}(3,:)>baseDims(3) ;
        roiCoords{iRoi}(:,outside_voxels) = [];
  end
  
  if ~isempty(roiCoords{iRoi})
    d.roiSize(iRoi) = size(roiCoords{iRoi},2);
    %remember which voxels are in which rois
    d.overlayRoiIndex{iRoi} = size(allRoisCoords,2)+(1:d.roiSize(iRoi));
    %put Roi coords together to simplify overlay computation
    allRoisCoords = [allRoisCoords roiCoords{iRoi}];

    %compute coordinates of all cube edges around each voxel
    edgesCoords = roiCoords{iRoi}([1 1 2 2 3 3],:);
    edgesCoords = reshape(repmat(edgesCoords,12,1),6,12*size(roiCoords{iRoi},2))';
    edgesCoords = edgesCoords+repmat(cubeEdgesCoords,size(roiCoords{iRoi},2),1);
    %find unique edges around  this roi
    [gridCoords, uniqueIndices] = unique(edgesCoords,'rows');
    doubleIndices = setdiff(1:size(edgesCoords,1),uniqueIndices);

    uniqueEdgeCoords = setdiff(edgesCoords,edgesCoords(doubleIndices,:),'rows');

    d.roiGridXcoords{iRoi} = gridCoords(:,1:2)';
    d.roiGridYcoords{iRoi} = gridCoords(:,3:4)';
    d.roiGridZcoords{iRoi} = gridCoords(:,5:6)';
    d.roiEdgeXcoords{iRoi} = uniqueEdgeCoords(:,1:2)';
    d.roiEdgeYcoords{iRoi} = uniqueEdgeCoords(:,3:4)';
    d.roiEdgeZcoords{iRoi} = uniqueEdgeCoords(:,5:6)';
    
    [d.roiGridXcoords{iRoi},d.roiGridYcoords{iRoi},d.roiGridZcoords{iRoi}] = convert2baseCoords(...
      d.roiGridXcoords{iRoi},d.roiGridYcoords{iRoi},d.roiGridZcoords{iRoi},d.base2scan);
    [d.roiEdgeXcoords{iRoi},d.roiEdgeYcoords{iRoi},d.roiEdgeZcoords{iRoi}] = convert2baseCoords(...
      d.roiEdgeXcoords{iRoi},d.roiEdgeYcoords{iRoi},d.roiEdgeZcoords{iRoi},d.base2scan);
    
    %Faces (opaque ROi)
    %coordinates of a cube around each voxel
    [d.roiFaceXcoords{iRoi}, d.roiFaceYcoords{iRoi}, d.roiFaceZcoords{iRoi}] = makeCubeFaces(roiCoords{iRoi},[],[],1,d.base2scan);
    %Faces (opaque ROi));

   else
      d.roiGridXcoords{iRoi} = [];
      d.roiGridYcoords{iRoi} = [];
      d.roiGridZcoords{iRoi} = [];
      d.roiEdgeXcoords{iRoi} = [];
      d.roiEdgeYcoords{iRoi} = [];
      d.roiEdgeZcoords{iRoi} = [];
      d.roiFaceXcoords{iRoi} = [];
      d.roiFaceYcoords{iRoi} = [];
      d.roiFaceZcoords{iRoi} = [];
      d.roiSize(iRoi) = 0;
      d.overlayRoiIndex{iRoi}=[];
   end
end
fprintf(1,'Done\n');

%find indices of unique voxels in all rois
allRoisCoords = unique(allRoisCoords','rows')';

%compute coordinates of a box in base space
allRoisBaseCoords = d.base2scan\[allRoisCoords;ones(1,size(allRoisCoords,2))];
minBaseCoords = floor(min(allRoisBaseCoords,[],2)');
maxBaseCoords = ceil(max(allRoisBaseCoords,[],2)');
if d.baseType==0
  minBaseCoords = max(minBaseCoords(1:3),[1 1 1]);    
  maxBaseCoords = min(maxBaseCoords(1:3),baseDims);   
end


%------------------------------------ Get Base data
baseColormap = gray(256);
baseClip = viewGet(thisView,'baseClip',baseNum);
baseGamma = viewGet(thisView,'baseGamma',baseNum);
switch(d.baseType)
  case 0
    for iDim=1:3
      otherDims = setdiff(1:3,iDim);
      for iSide=1:2
        switch iSide
          case 1
            sliceNum = minBaseCoords(iDim);
          case 2
            sliceNum = maxBaseCoords(iDim);
        end
        [baseImage,baseCoords{iDim,iSide},baseCoordsHomogeneous{iDim,iSide}] = ...
          getBaseSlice(thisView,sliceNum,iDim,0,baseNum,d.baseType);
        d.baseRGB{iDim,iSide} = rescale2rgb(baseImage,baseColormap,baseClip,baseGamma);
        thisMinCoords = minBaseCoords(otherDims);
        thisMaxCoords = maxBaseCoords(otherDims);
        d.baseRGB{iDim,iSide} = d.baseRGB{iDim,iSide}(thisMinCoords(1):thisMaxCoords(1),thisMinCoords(2):thisMaxCoords(2),:);
        d.baseRGB{iDim,iSide}(:,end+1,:) = d.baseRGB{iDim,iSide}(:,end,:); %Add one column
        d.baseRGB{iDim,iSide}(end+1,:,:) = d.baseRGB{iDim,iSide}(end,:,:); %Add one row
        baseCoords{iDim,iSide} = baseCoords{iDim,iSide}(thisMinCoords(1):thisMaxCoords(1),thisMinCoords(2):thisMaxCoords(2),:);
        baseCoords{iDim,iSide}(:,end+1,:) = baseCoords{iDim,iSide}(:,end,:); %Add one column
        baseCoords{iDim,iSide}(:,end,otherDims(2)) = baseCoords{iDim,iSide}(:,end,otherDims(2))+1; %Add one column
        baseCoords{iDim,iSide}(end+1,:,:) = baseCoords{iDim,iSide}(end,:,:); %Add one row
        baseCoords{iDim,iSide}(end,:,otherDims(1)) = baseCoords{iDim,iSide}(end,:,otherDims(1))+1; %Add one row
        baseCoords{iDim,iSide}(:,:,otherDims) = baseCoords{iDim,iSide}(:,:,otherDims) - .5; %shift the coords by .5 so that each square will be centered on the integer coordinate
        switch iSide
          case 1
            baseCoords{iDim,iSide}(:,:,iDim) = baseCoords{iDim,iSide}(:,:,iDim) - .5; %on the front side, shift the slice backwards
          case 2
            baseCoords{iDim,iSide}(:,:,iDim) = baseCoords{iDim,iSide}(:,:,iDim) + .5; %and on the back side, shift it forwards
        end
      end
    end
    
  case 1
    baseImage = baseImage(minBaseCoords(1):maxBaseCoords(1),minBaseCoords(2):maxBaseCoords(2));
    d.baseRGB{3,1} = rescale2rgb(baseImage,baseColormap,baseClip,baseGamma);
    d.baseRGB{3,2} = d.baseRGB{3,1};
    [baseCoords{3,1}(:,:,1),baseCoords{3,1}(:,:,2),baseCoords{3,1}(:,:,3)] = ndgrid(minBaseCoords(1):maxBaseCoords(1),minBaseCoords(2):maxBaseCoords(2),minBaseCoords(3)-.5);
    [baseCoords{3,2}(:,:,1),baseCoords{3,2}(:,:,2),baseCoords{3,2}(:,:,3)] = ndgrid(minBaseCoords(1):maxBaseCoords(1),minBaseCoords(2):maxBaseCoords(2),maxBaseCoords(3)+.5);
    
  case 2
    %get the surface coordinates
    baseCoordMap = viewGet(thisView,'baseCoordMap',baseNum);
    d.surfFaces = baseCoordMap.tris;
    d.innerVertices = permute(baseCoordMap.innerCoords,[2 4 1 3]);
    d.outerVertices = permute(baseCoordMap.outerCoords,[2 4 1 3]);
    corticalDepths = str2num(num2str(baseCoordMap.corticalDepths));
    clear('baseCoordMap');
    d.corticalDepth = str2num(num2str(viewGet(thisView,'corticalDepth')));
    baseData = viewGet(thisView,'baseData',baseNum);
% % %     %if we get the surface coordinates using viewGet(thisView,'surfaceData'), we need to
% % %     %swap X and Y (Not sure why I have to do that and how general it is, but it works on the present data)
% % %     tempVtcs = d.innerVertices(:,1);
% % %     d.innerVertices(:,1) = d.innerVertices(:,2);
% % %     d.innerVertices(:,2) = tempVtcs;
% % %     tempVtcs = d.outerVertices(:,1);
% % %     d.outerVertices(:,1) = d.outerVertices(:,2);
% % %     d.outerVertices(:,2) = tempVtcs;
    %find vtcs that are in the box (at mid cortical depth)
    midVertices = (d.outerVertices+d.innerVertices)/2;
    verticesInBox = midVertices(:,1)>=minBaseCoords(1) & midVertices(:,1)<=maxBaseCoords(1) &...
                    midVertices(:,2)>=minBaseCoords(2) & midVertices(:,2)<=maxBaseCoords(2) &...
                    midVertices(:,3)>=minBaseCoords(3) & midVertices(:,3)<=maxBaseCoords(3);
    %remove vtcs outside ROI
    d.innerVertices = d.innerVertices(verticesInBox,:);
    d.outerVertices = d.outerVertices(verticesInBox,:);
    baseData = baseData(verticesInBox);
    %replace ones by its 'vertex in the box' number
    verticesInBox = double(verticesInBox);     %need to convert to double first     
    verticesInBox(verticesInBox>0)=double(1:nnz(verticesInBox));
    %remove any face that involves vtcs outside the box
    d.surfFaces = verticesInBox(d.surfFaces);
    d.surfFaces = d.surfFaces(all(d.surfFaces>0,2),:);
    d.baseRGB = rescale2rgb(baseData,baseColormap,baseClip,baseGamma);
    d.baseRGB = permute(d.baseRGB,[2 3 1]);
end

%--------------------------------Mask and Alpha overlay data -------------------------------
if ~isempty(currentAnalysis)
   
  d.nOverlays = length(overlayList);
  d.alphaOverlayList = zeros(1,d.nOverlays);
  for iOverlay = 1:d.nOverlays
    d.overlayClip(iOverlay,:) = str2num(num2str(viewGet(thisView,'overlayClip',overlayList(iOverlay)))); 
    thisNum = viewGet(thisView,'overlayNum',viewGet(thisView,'alphaOverlay',overlayList(iOverlay)));
    if ~isempty(thisNum)
      d.alphaOverlayList(iOverlay) = thisNum;
      d.alphaOverlayClip(iOverlay,:) = str2num(num2str(viewGet(thisView,'overlayClip',d.alphaOverlayList(iOverlay))));
    else
      d.alphaOverlayClip(iOverlay,:) = [-1 1];
    end
  end
  d.overlayCurrentClip = d.overlayClip;
  d.alphaOverlayCurrentClip = d.alphaOverlayClip;

  %First, get a mask of non-zero voxel representing the current overlay display
  %This is taken from computeOverlay.m
  fprintf(1,'Computing overlays...');
  [d.overlayMaskData, d.overlayData]= maskOverlay(thisView,[overlayList d.alphaOverlayList],scanNum);
  d.overlayMaskData = d.overlayMaskData{1};
  d.overlayData = d.overlayData{1};
  d.overlayMaskData(isnan(d.overlayMaskData))=0;
  if d.baseType == 1 
    d.overlayData = getBaseSpaceOverlay(thisView,d.overlayData,scanNum,baseNum,[], corticalDepthBins, rotateAngle);
    d.overlayMaskData = getBaseSpaceOverlay(thisView,double(d.overlayMaskData),scanNum,baseNum,'nearest', corticalDepthBins, rotateAngle);
    d.overlayData = d.overlayData(:,:,corticalDepthIndices,:);
    d.overlayMaskData = d.overlayMaskData(:,:,corticalDepthIndices,:);
    d.overlayMaskData(isnan(d.overlayMaskData))=0;
  end

  d.overlayMaskData = logical(d.overlayMaskData);

  fprintf(1,'Done\n');
   %make alpha map
   if ~any(d.alphaOverlayList)
      d.overlayAlpha = ones(size(d.overlayMaskData));
      d.overlayAlphaMaskData = ones(size(d.overlayMaskData));
      d.overlayAlphaData = zeros(size(d.overlayMaskData));
   else
      fprintf(1,'Computing alpha overlay...');
      
      d.overlayAlphaData = d.overlayData(:,:,:,d.nOverlays+1:end);
      d.overlayData = d.overlayData(:,:,:,1:d.nOverlays);
      
      d.overlayAlphaMaskData = d.overlayMaskData(:,:,:,d.nOverlays+1:end);
      d.overlayMaskData = d.overlayMaskData(:,:,:,1:d.nOverlays);
 
      % get the range of the alpha overlay 
      cOverlay=0;
      for iOverlay = d.alphaOverlayList;
        cOverlay=cOverlay+1;
        if ~d.alphaOverlayList(cOverlay)
          d.overlayAlpha(:,:,:,cOverlay) = ones(overlayDims);
          d.overlayAlphaData(:,:,:,cOverlay) = zeros(overlayDims);
          d.overlayAlphaMaskData(:,:,:,cOverlay) = ones(overlayDims);
        else
          range = viewGet(thisView,'overlayClip',iOverlay);
          % handle setRangeToMax (to debug)
          if strcmp(viewGet(thisView,'overlayCtype',d.alphaOverlayList),'setRangeToMax')
            clip = viewGet(thisView,'overlayClip',iOverlay);
           maxRange = max(clip(1),min(d.overlayAlphaData(d.overlayAlphaMaskData)));
           if ~isempty(maxRange),range(1) = maxRange;end
           minRange = min(max(d.overlayAlphaData(d.overlayAlphaMaskData)),clip(2));
           if ~isempty(minRange),range(2) = minRange;end
          end
          % now compute the alphaOverlay as a number from
          % 0 to 1 of the range
          thisAlphaOverlay = (d.overlayAlphaData(:,:,:,cOverlay)-range(1))./diff(range);
%           thisAlphaOverlay(thisAlphaOverlay>1) = 1;
%           thisAlphaOverlay(thisAlphaOverlay<0) = 0;
          alphaOverlayExponent = viewGet(thisView,'alphaOverlayExponent');
          if alphaOverlayExponent<0   %JB if the alpha overlay exponent is negative, set it positive and take the complementary values to 1 for the alpha map
            alphaOverlayExponent = -alphaOverlayExponent;
            thisAlphaOverlay = 1-thisAlphaOverlay;
          end
          d.overlayAlpha(:,:,:,cOverlay) = thisAlphaOverlay.^alphaOverlayExponent;
        end
      end
      d.overlayAlpha(isnan(d.overlayAlpha)) = 0;
   end   
   fprintf(1,'Done\n');

end


%------------------------------------ Reduce Overlay Data  ------------------------------------------------%

if ~isempty(currentAnalysis)
  %find indices of all scan voxels within the base space box
  [scanXCoords,scanYCoords,scanZCoords] = ndgrid(1:overlayDims(1),1:overlayDims(2),1:overlayDims(3));
  baseScanCoords = d.base2scan\[reshape(cat(4,scanXCoords,scanYCoords,scanZCoords),prod(overlayDims),3)';ones(1,prod(overlayDims))];
  scanVoxelsInBox = baseScanCoords(1,:)>=minBaseCoords(1) & baseScanCoords(1,:)<=maxBaseCoords(1) & ...
                    baseScanCoords(2,:)>=minBaseCoords(2) & baseScanCoords(2,:)<=maxBaseCoords(2) & ...
                    baseScanCoords(3,:)>=minBaseCoords(3) & baseScanCoords(3,:)<=maxBaseCoords(3);

  %restrict data to these voxels
  d.overlayAlphaData = permute(d.overlayAlphaData,[4 1 2 3]);
  d.overlayAlphaData = d.overlayAlphaData(:,scanVoxelsInBox)';
  d.overlayAlpha = permute(d.overlayAlpha,[4 1 2 3]);
  d.overlayAlpha = d.overlayAlpha(:,scanVoxelsInBox)';
  d.overlayData = permute(d.overlayData,[4 1 2 3]);
  d.overlayData = d.overlayData(:,scanVoxelsInBox)';
  d.overlayMaskData = permute(d.overlayMaskData,[4 1 2 3]);
  d.overlayMaskData = d.overlayMaskData(:,scanVoxelsInBox)';
  d.overlayAlphaMaskData = permute(d.overlayAlphaMaskData,[4 1 2 3]);
  d.overlayAlphaMaskData = d.overlayAlphaMaskData(:,scanVoxelsInBox)';
  %keep scan coordinates  of these voxels
  d.overlayScanCoords = [scanXCoords(scanVoxelsInBox);  scanYCoords(scanVoxelsInBox); scanZCoords(scanVoxelsInBox)];  

  %Now remove voxels that will never be displayed, 
  %either because they're masked in all overlays
  voxelsToKeep = any(d.overlayMaskData,2);
  %because they're masked in all alpha overlays
  voxelsToKeep = voxelsToKeep & any(d.overlayAlphaMaskData,2);
  %because they're zero in all alpha overlays
  voxelsToKeep = voxelsToKeep & any(d.overlayAlpha,2);
  %or because at least they're NaN in at least one overlay
  %(this condition could be made less restrictive and be changed to: they're NaN in all overlays)
  voxelsToKeep = voxelsToKeep & ~all(isnan(d.overlayData),2);
  d.overlayAlphaData = d.overlayAlphaData(voxelsToKeep,:);
  d.overlayAlpha = d.overlayAlpha(voxelsToKeep,:);
  d.overlayData = d.overlayData(voxelsToKeep,:);
  d.overlayMaskData = d.overlayMaskData(voxelsToKeep,:);
  d.overlayAlphaMaskData = d.overlayAlphaMaskData(voxelsToKeep,:);
  d.overlayScanCoords = d.overlayScanCoords(:,voxelsToKeep);
  
  %make a copy of the masks
  d.overlayClippedMask = d.overlayMaskData;
  d.overlayAlphaClippedMask = d.overlayAlphaMaskData;

  %remember which data points are in which ROI
  for iRoi = 1:nRois
    [~,d.overlayRoiIndex{iRoi}] = intersect(d.overlayScanCoords',roiCoords{iRoi}','rows');
  end
  [~,d.overlayAllRoisIndex] = intersect(d.overlayScanCoords',allRoisCoords','rows');

  cOverlay = 0;
  d.overlayRGB = NaN([size(d.overlayData,1) 1 3 size(d.overlayData,2)]);
  for iOverlay=overlayList
    cOverlay = cOverlay+1;
    %d.colorMap = viewGet(thisView,'colormap',iOverlay);
    colorMap = analysis.overlays(iOverlay).colormap;
    colorRange = str2num(num2str(viewGet(thisView,'overlayClip',iOverlay))); %str2num(num2str()) is because the edit box rounds the values  
                                                                            %and causes problems if the min/max of the sliders is not rounded
    d.overlayRGB(:,:,:,cOverlay) = rescale2rgb(d.overlayData(:,cOverlay),colorMap,colorRange);
  end
  d.overlayRGB = permute(d.overlayRGB,[1 3 4 2]);
  
      
end
      
toc

%---------------------- construct figure -----------------%
hFigure = selectGraphWin;
% turn off menu/title etc.
set(hFigure,'NumberTitle','off','Name','RenderRois3D','Toolbar','Figure','color',[.7 .7 .7],'Renderer','OpenGL');
roisButtonsBottom = .5;
overlayButtonsBottom = .04;
baseButtonsBottom = .8;
overlaySlidersLeft=.02;
rotateViewLeft = .6;
margin=.02;
sliderWidth = .1;
stringWidth = .03;

d.hMainAxis = gca;
set(d.hMainAxis,'Alim',[0 1],'color','none')


%Rois plot and control buttons
hold on
d.roiIsEmpty = [false(1,length(roiName)) true];
for iRoi = 1:nRois
   temp = patch(d.roiEdgeXcoords{iRoi},d.roiEdgeYcoords{iRoi},d.roiEdgeZcoords{iRoi},roiColor(iRoi,:),'edgecolor',roiColor(iRoi,:),'LineWidth',2,'parent',d.hMainAxis);
   if ~isempty(temp)
      d.hRoi(iRoi) = temp;
   else
     d.roiIsEmpty(iRoi)=true;
   end
end
uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','popupmenu',       'Position',[margin roisButtonsBottom+.2 .1 .03],...
      'Value',1, 'String', {'Perimeter','Grid','Opaque'},'Callback',{@drawRois,d});
d.hRoiList = uicontrol('Parent',hFigure, 'Unit','normalized',               'Position',[margin roisButtonsBottom .15 .2],...
  'Style','listBox', 'min',0,'max',length(roiName),'Value',1:length(roiName), 'String', roiName,'Callback',@showRois);
d.hShowROIs = uicontrol('Parent',hFigure, 'Unit','normalized',         'Position',[margin roisButtonsBottom+.23 .15 .03],...
         'Style','Checkbox','Callback',@showRois, 'String','Show ROIs','value',1);

%Overlay control buttons
if ~isempty(currentAnalysis)
  overlayNames = viewGet(thisView,'overlayNames');
  overlayNames = overlayNames(overlayList);
   set(hFigure,'Name',['RenderRois3D - ' viewGet(thisView,'homeDir') ' - ' viewGet(thisView,'analysisName')]);

   overlayModes = {'Show all voxels','Restrict to all ROIs','Restrict to selected ROIs'};
   d.hOverlayRoiMode = uicontrol('Parent',hFigure, 'Unit','normalized',     'Position',[margin overlayButtonsBottom+.36 .15 .03],...
      'Style','popupmenu', 'Value',1, 'String', overlayModes,'Callback',@recomputeOverlay);
    
   d.hOverlayList = uicontrol('Parent',hFigure, 'Unit','normalized',        'Position',[margin overlayButtonsBottom+.16 .15 .2],...
      'Style','listBox', 'min',0,'max',length(overlayNames),'Value',1:length(overlayNames), 'String', overlayNames,'Callback',@recomputeOverlay);
   
   %Overlay clip buttons
   uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','text',       'Position',[margin overlayButtonsBottom+.12 stringWidth .03],'String','Min.');
   d.hOverlayEditMin = uicontrol('Parent',hFigure, 'Unit','normalized',   'Position',[overlaySlidersLeft+stringWidth+sliderWidth overlayButtonsBottom+.13 .03 .03],...
      'Style','edit', 'String',num2str(d.overlayClip(1)),'Callback',@recomputeOverlay);
   d.hOverlaySliderMin = uicontrol('Parent',hFigure, 'Unit','normalized', 'Position',[overlaySlidersLeft+stringWidth overlayButtonsBottom+.12 sliderWidth .03],...
      'Style','slider', 'Min', d.overlayClip(1,1),'Max',d.overlayClip(1,2), 'value',d.overlayClip(1,1),'Callback',@recomputeOverlay);

   uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','text',       'Position',[margin overlayButtonsBottom+.1 stringWidth .03],'String','Max.');
   d.hOverlayEditMax = uicontrol('Parent',hFigure, 'Unit','normalized',   'Position',[overlaySlidersLeft+stringWidth+sliderWidth overlayButtonsBottom+.11 .03 .03],...
      'Style','edit', 'String',num2str(d.overlayClip(1,2)),'Callback',@recomputeOverlay);
   d.hOverlaySliderMax = uicontrol('Parent',hFigure, 'Unit','normalized', 'Position',[overlaySlidersLeft+stringWidth overlayButtonsBottom+.10 sliderWidth .03],...
      'Style','slider', 'Min', d.overlayClip(1,1),'Max',d.overlayClip(1,2), 'value',d.overlayClip(1,2),'Callback',@recomputeOverlay);
   if d.nOverlays>1
    set([d.hOverlayEditMin d.hOverlaySliderMin d.hOverlayEditMax d.hOverlaySliderMax],'enable','off')
   end
   
   %Overlay alpha buttons
   uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','text',         'Position',[margin overlayButtonsBottom+.08 stringWidth .03],'String','Alpha');
   d.hOverlayAlphaEdit = uicontrol('Parent',hFigure, 'Unit','normalized',   'Position',[overlaySlidersLeft+stringWidth+sliderWidth overlayButtonsBottom+.09 .03 .03],...
      'Style','edit', 'String',overlayAlpha,'Callback',@recomputeOverlay);
   d.hOverlayAlphaSlider = uicontrol('Parent',hFigure, 'Unit','normalized', 'Position',[overlaySlidersLeft+stringWidth overlayButtonsBottom+.08 sliderWidth .03],...
      'Style','slider', 'Min', 0,'Max',1, 'value',overlayAlpha,'Callback',@recomputeOverlay);
    
   %Alpha overlay buttons
   if any(d.alphaOverlayList)
     d.hOverlayAlphaCheck = uicontrol('Parent',hFigure, 'Unit','normalized',    'Position',[margin overlayButtonsBottom+.05 .15 .03],...
      'Style','checkbox', 'Min', 0,'Max',1, 'Value',1, 'String','Use Alpha Overlay','Callback',@recomputeOverlay);
     
     uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','text',           'Position',[margin overlayButtonsBottom+.02 stringWidth .03],'String','Min.');
     d.hAlphaOverlayEditMin = uicontrol('Parent',hFigure, 'Unit','normalized',  'Position',[overlaySlidersLeft+stringWidth+sliderWidth overlayButtonsBottom+.03 .03 .03],...
        'Style','edit', 'String',num2str(d.alphaOverlayClip(1,1)),'Callback',@recomputeOverlay);
     d.hAlphaOverlaySliderMin = uicontrol('Parent',hFigure, 'Unit','normalized', 'Position',[overlaySlidersLeft+stringWidth overlayButtonsBottom+.02 sliderWidth .03],...
        'Style','slider', 'Min', d.alphaOverlayClip(1,1),'Max',d.alphaOverlayClip(1,2), 'value',d.alphaOverlayClip(1,1),'Callback',@recomputeOverlay);
     
     uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','text',           'Position',[margin overlayButtonsBottom+.00 stringWidth .03],'String','Max.');
     d.hAlphaOverlayEditMax = uicontrol('Parent',hFigure, 'Unit','normalized',  'Position',[overlaySlidersLeft+stringWidth+sliderWidth overlayButtonsBottom+.01 .03 .03],...
        'Style','edit', 'String',num2str(d.alphaOverlayClip(1,2)),'Callback',@recomputeOverlay);
     d.hAlphaOverlaySliderMax = uicontrol('Parent',hFigure, 'Unit','normalized', 'Position',[overlaySlidersLeft+stringWidth overlayButtonsBottom+.00 sliderWidth .03],...
        'Style','slider', 'Min', d.alphaOverlayClip(1,1),'Max',d.alphaOverlayClip(1,2), 'value',d.alphaOverlayClip(1,2),'Callback',@recomputeOverlay);
     if d.nOverlays>1
      set([d.hAlphaOverlayEditMin d.hAlphaOverlaySliderMin d.hAlphaOverlayEditMax d.hAlphaOverlaySliderMax],'enable','off')
     end
   end
   
   set(hFigure,'userdata',d);
   recomputeOverlay(d.hOverlayAlphaEdit);
   
   %this has to be after calling recomputeOverlay, because we want to use the hande of the overlay cubes
   d = get(hFigure,'userdata');
   set(d.hOverlay,'ambientStrength',.5);  %not too bright, not too dark
   set(d.hOverlay,'specularStrength',.3); %little reflectance
   uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','Checkbox', 'Position',[margin overlayButtonsBottom+.39 .15 .03],...
         'Callback',{@makeVisible,d.hOverlay}, 'String','Show Overlay','value',1);

end

% plot Base slices
switch(d.baseType)
  case {0,1}
    for iDim = 1:3
      for iSide = 1:2
        if ~isempty(baseCoords{iDim,iSide})
          d.hSurface(iDim,iSide) = surface(baseCoords{iDim,iSide}(:,:,1),baseCoords{iDim,iSide}(:,:,2),baseCoords{iDim,iSide}(:,:,3),d.baseRGB{iDim,iSide}, 'EdgeColor','none','parent',d.hMainAxis);
%         else
%           d.hSurface(iDim,iSide) = [];
        end
      end
    end
    uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','checkbox',   'Position',[margin baseButtonsBottom+.1 .15 .03],...
         'Min', 0,'Max',1, 'Value',1, 'String', 'Show Base', 'callback',{@makeBaseVisible,d.hSurface});
  case 2
    sliderStep(1) = corticalDepths(2) - corticalDepths(1);
    sliderStep(2) = 3*sliderStep(1);
    %depth control
    uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','text',       'Position',[margin baseButtonsBottom+.09 stringWidth .03],'String','Depth');
    uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','text',       'Position',[margin baseButtonsBottom+.02 stringWidth .03],'String','Depth');
    d.hInnerSurfEdit = uicontrol('Parent',hFigure, 'Unit','normalized',    'Position',[margin+stringWidth+sliderWidth baseButtonsBottom+.1 .02 .03],...
      'Style','edit', 'String',num2str(d.corticalDepth(1)),'Callback',@recomputeSurface);
    d.hOuterSurfEdit = uicontrol('Parent',hFigure, 'Unit','normalized',    'Position',[margin+stringWidth+sliderWidth baseButtonsBottom+.03 .02 .03],...
      'Style','edit', 'String',num2str(d.corticalDepth(2)),'Callback',@recomputeSurface);
    d.hInnerSurfSlider = uicontrol('Parent',hFigure, 'Unit','normalized',  'Position',[margin+stringWidth baseButtonsBottom+.09 sliderWidth .03],...
      'Style','slider', 'Min', corticalDepths(1),'Max',corticalDepths(end), 'value',d.corticalDepth(1),'sliderStep',sliderStep,'Callback',@recomputeSurface);
    d.hOuterSurfSlider = uicontrol('Parent',hFigure, 'Unit','normalized',  'Position',[margin+stringWidth baseButtonsBottom+.02 sliderWidth .03],...
      'Style','slider', 'Min', corticalDepths(1),'Max',corticalDepths(end), 'value',d.corticalDepth(2),'sliderStep',sliderStep,'Callback',@recomputeSurface);
    
    %transparency control
    uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','text',       'Position',[margin baseButtonsBottom+.07 stringWidth .03],'String','Alpha');
    uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','text',       'Position',[margin baseButtonsBottom stringWidth .03],'String','Alpha');
    d.hInnerAlphaEdit = uicontrol('Parent',hFigure, 'Unit','normalized',   'Position',[margin+stringWidth+sliderWidth baseButtonsBottom+.08 .02 .03],...
      'Style','edit', 'String',num2str(1),'Callback',@recomputeSurface);
    d.hOuterAlphaEdit = uicontrol('Parent',hFigure, 'Unit','normalized',   'Position',[margin+stringWidth+sliderWidth baseButtonsBottom+.01 .02 .03],...
      'Style','edit', 'String',num2str(1),'Callback',@recomputeSurface);
    d.hInnerAlphaSlider = uicontrol('Parent',hFigure, 'Unit','normalized', 'Position',[margin+stringWidth baseButtonsBottom+.07 sliderWidth .03],...
      'Style','slider', 'Min', 0,'Max',1, 'value',1,'sliderStep',[.1 .3],'Callback',@recomputeSurface);
    d.hOuterAlphaSlider = uicontrol('Parent',hFigure, 'Unit','normalized', 'Position',[margin+stringWidth baseButtonsBottom sliderWidth .03],...
      'Style','slider', 'Min', 0,'Max',1, 'value',1,'sliderStep',[.1 .3],'Callback',@recomputeSurface);
   
    set(hFigure,'userdata',d);
    recomputeSurface(d.hOuterSurfSlider)
    recomputeSurface(d.hInnerSurfSlider)
    d = get(hFigure,'userdata');
%    shading(d.hMainAxis,'flat')
    uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','checkbox','Position',[margin baseButtonsBottom+.05 .15 .03],...
         'Min', 0,'Max',1, 'Value',1, 'String', 'Show Outer Surface', 'callback',{@makeVisible,d.hOuterSurface});
    uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','checkbox','Position',[margin baseButtonsBottom+.12 .15 .03],...
         'Min', 0,'Max',1, 'Value',1, 'String', 'Show Inner Surface', 'callback',{@makeVisible,d.hInnerSurface});
    
end
if ismember(d.baseType,[0 2])
end

daspect(1./baseVoxelSize)
view(3);
axis tight vis3d;
lighting phong
grid on;
% if ~isempty(currentAnalysis)
%    colormap(d.colorMap);
%    caxis(d.colorRange);
%    colorbar('Units','normalized', 'OuterPosition',[.91 .5 .08 .45 ]);
% end

xlabel('X');
ylabel('Y');
zlabel('Z');
set(gca,'Xdir','reverse'); %the X axis goes from right to left(?)
d.reversedX = true;
uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','checkbox', 'Position',[rotateViewLeft+.1 .04 .15 .03],...
         'Min', 0,'Max',1, 'Value',d.reversedX, 'String', 'Reverse X Axis','Callback',{@reverseX});

%rotation control
hRotate = rotate3d;
set(hRotate,'Enable','on');
[azimuth,elevation] = view;
rotateViewSliderLength = .90-rotateViewLeft;
uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','text', 'Position',[rotateViewLeft .03 .05 .03],'String', 'Azimuth');
d.hAzimuthEdit = uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','edit', 'Position',[rotateViewLeft .01 .05 .03],...
      'String',num2str(azimuth));
d.hAzimuthSlider = uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','slider', 'Position',[rotateViewLeft+.05 .01 rotateViewSliderLength .03],...
      'Min', -360,'Max',360, 'value',azimuth);
uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','text', 'Position',[.94 .07+rotateViewSliderLength .05 .03],'String', 'Elevation');
d.hElevationEdit = uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','edit', 'Position',[.94 .05+rotateViewSliderLength .05 .03 ],...
      'String',num2str(elevation));
d.hElevationSlider = uicontrol('Parent',hFigure, 'Unit','normalized', 'Style','slider', 'Position',[.95 .05 .03 rotateViewSliderLength],...
      'Min', -90,'Max',90, 'value',elevation);
set(d.hAzimuthSlider,'Callback',@changeView);
set(d.hElevationSlider,'Callback',@changeView);
set(d.hAzimuthEdit,'Callback',@changeView);
set(d.hElevationEdit,'Callback',@changeView);
set(hRotate,'ActionPostCallback',@changeView);

%set 2 lights at opposite sides of the object
h_light1 = camlight(90,45);
h_light2 = camlight(-90,-45);

set(hFigure,'userdata',d);
[azimuth,elevation] = view(d.hMainAxis);
setBaseVisibility(d,azimuth,elevation);


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%---------------------------------------------------- change and/or diplays the elevation and azimuth after mouse rotate (rotate3d)
function changeView(handle,eventdata)

if isnumeric(handle) %if the function has been called from one of the figure uicontrols
  hFigure= get(handle,'parent');
  d = get(hFigure,'userdata');
  azimuth = round(get(d.hAzimuthSlider,'Value'));
  elevation = round(get(d.hElevationSlider,'Value'));
  switch(handle)
    case d.hAzimuthSlider
       azimuth = round(get(handle,'Value'));
    case d.hElevationSlider
       elevation = round(get(handle,'Value'));
    case d.hAzimuthEdit
        azimuth = round(str2num(get(handle,'String')));
    case d.hElevationEdit
       elevation= round(str2num(get(handle,'String')));
  end
  view(azimuth,elevation);
else %if the function has been called after rotating the axes with the mouse
  hFigure= get(eventdata.Axes,'parent');
  d = get(hFigure,'userdata');
  [azimuth,elevation] = view(d.hMainAxis);
  azimuth = round(azimuth);
  elevation = round(elevation);
end

set(d.hAzimuthSlider,'Value',azimuth);
set(d.hElevationSlider,'Value',elevation);
set(d.hAzimuthEdit,'String',num2str(azimuth));
set(d.hElevationEdit,'String',num2str(elevation));

%show or hide base slices depending on the orientation of the view
setBaseVisibility(d,azimuth,elevation);

return

%----------------------------------------------------  sets Base visibility depending on view
function setBaseVisibility(d,azimuth,elevation)

if ismember(d.baseType,[0 1])
  if d.reversedX
    xNegative = -sin(azimuth*pi/180); %this is negative if the viewpoint is X-positive
  else
    xNegative = sin(azimuth*pi/180); %this is positive if the viewpoint is X-positive
  end
  if xNegative>0 
    set(d.hSurface(1,1),'visible','on');
    set(d.hSurface(1,2),'visible','off');
  else
    set(d.hSurface(1,1),'visible','off');
    set(d.hSurface(1,2),'visible','on');
  end    
  if -cos(azimuth*pi/180)>0 %this is positive if the viewpoint is Y-positive
    set(d.hSurface(2,1),'visible','on');
    set(d.hSurface(2,2),'visible','off');
  else
    set(d.hSurface(2,1),'visible','off');
    set(d.hSurface(2,2),'visible','on');
  end    
  if sin(elevation*pi/180)>0 %this is positive if the viewpoint is Z-positive
    set(d.hSurface(3,1),'visible','on');
    set(d.hSurface(3,2),'visible','off');
  else
    set(d.hSurface(3,1),'visible','off');
    set(d.hSurface(3,2),'visible','on');
  end    
end

return;


%---------------------------------------------------- reverse X axis
function reverseX(handle,eventdata)
hFigure= get(handle,'parent');
d = get(hFigure,'userdata');

d.reversedX = get(handle,'value');
if d.reversedX
   set(gca,'xdir','reverse');
else
   set(gca,'xdir','normal');
end
if ~d.baseType
%switch which coronal base slice is displayed 
  if strcmp(get(d.hSurface(1,1),'visible'),'on')
    set(d.hSurface(1,1),'visible','off');
    set(d.hSurface(1,2),'visible','on');
  else
    set(d.hSurface(1,1),'visible','on');
    set(d.hSurface(1,2),'visible','off');
  end
end

set(hFigure,'userdata',d);
return;


%---------------------------------------------------- makes objects visible or invisible
function makeVisible(handle,eventdata,hObject)

if get(handle,'value')
   set(hObject,'visible','on');
else
   set(hObject,'visible','off');
end

return;

%---------------------------------------------------- makes objects visible or invisible
function makeBaseVisible(handle,eventdata,hObject)

if get(handle,'value')
  hFigure= get(handle,'parent');
  d = get(hFigure,'userdata');
  [azimuth,elevation] = view(d.hMainAxis);
  setBaseVisibility(d,azimuth,elevation);
else
  set(hObject,'visible','off');
end

return;

%---------------------------------------------------- draws Rois
function drawRois(handle,eventdata,d)

hFigure= get(handle,'parent');
d = get(hFigure,'userdata');

list = get(handle,'String');
switch(list{get(handle,'Value')})
   case 'Perimeter'
      Xcoords = d.roiEdgeXcoords;
      Ycoords = d.roiEdgeYcoords;
      Zcoords = d.roiEdgeZcoords;
   case 'Grid'
      Xcoords = d.roiGridXcoords;
      Ycoords = d.roiGridYcoords;
      Zcoords = d.roiGridZcoords;
   case 'Opaque'
      Xcoords = d.roiFaceXcoords;
      Ycoords = d.roiFaceYcoords;
      Zcoords = d.roiFaceZcoords;
end

for iRoi = 1:length(Xcoords)
   if ~isempty(Xcoords{iRoi})
      set(d.hRoi(iRoi),'Xdata',Xcoords{iRoi},'Ydata',Ycoords{iRoi},'Zdata',Zcoords{iRoi});
   end
end

%---------------------------------------------------- show Rois
function showRois(handle,eventdata)

hFigure= get(handle,'parent');
d = get(hFigure,'userdata');

roisToShow = false(1,length(d.roiIsEmpty));
roisToShow(get(d.hRoiList,'value'))=true;

if get(d.hShowROIs,'value')
  set(d.hRoi(~d.roiIsEmpty & roisToShow),'visible','on')
  set(d.hRoi(~d.roiIsEmpty & ~roisToShow),'visible','off')
else
  set(d.hRoi(~d.roiIsEmpty),'visible','off')
end

if isfield(d,'hOverlayRoiMode')
  switch(handle)
    case d.hRoiList
      overlayRoiModes = get(d.hOverlayRoiMode,'String');
      overlayRoiMode = overlayRoiModes{get(d.hOverlayRoiMode,'Value')};
      if strcmp(overlayRoiMode,'Restrict to selected ROIs')
        recomputeOverlay(handle);
      end
  end
end
%---------------------------------------------------- compute Surface coordinates
function recomputeSurface(handle,eventdata)

hFigure= get(handle,'parent');
d = get(hFigure,'userdata');

switch(handle)
  case d.hInnerSurfEdit
    thisCorticalDepth = str2num(get(handle,'string'));
    set(d.hInnerSurfSlider,'value',thisCorticalDepth);
  case d.hOuterSurfEdit
    thisCorticalDepth = str2num(get(handle,'string'));
    set(d.hOuterSurfSlider,'value',thisCorticalDepth);
  case d.hInnerSurfSlider
    thisCorticalDepth = get(handle,'value');
    set(d.hInnerSurfEdit,'string',num2str(thisCorticalDepth));
  case d.hOuterSurfSlider
    thisCorticalDepth = get(handle,'value');
    set(d.hOuterSurfEdit,'string',num2str(thisCorticalDepth));
  case d.hInnerAlphaEdit
    thisAlpha = str2num(get(handle,'string'));
    set(d.hInnerAlphaSlider,'value',thisAlpha);
  case d.hOuterAlphaEdit
    thisAlpha = str2num(get(handle,'string'));
    set(d.hOuterAlphaSlider,'value',thisAlpha);
  case d.hInnerAlphaSlider
    thisAlpha = get(handle,'value');
    set(d.hInnerAlphaEdit,'string',num2str(thisAlpha));
  case d.hOuterAlphaSlider
    thisAlpha = get(handle,'value');
    set(d.hOuterAlphaEdit,'string',num2str(thisAlpha));
end
    
switch(handle)
  
  case {d.hInnerSurfEdit, d.hInnerSurfSlider}
    vertices = d.innerVertices + thisCorticalDepth*(d.outerVertices-d.innerVertices);
    if isfield(d,'hInnerSurface')
      set(d.hInnerSurface,'vertices',vertices);
    else
      d.hInnerSurface = patch('vertices', vertices, 'faces', d.surfFaces,'FaceVertexCData', d.baseRGB,'facecolor','interp','edgecolor','none','parent',d.hMainAxis);
    end
    
  case {d.hOuterSurfEdit, d.hOuterSurfSlider}
    vertices = d.innerVertices + thisCorticalDepth*(d.outerVertices-d.innerVertices);
    if isfield(d,'hOuterSurface')
      set(d.hOuterSurface,'vertices',vertices);
    else
      d.hOuterSurface = patch('vertices', vertices, 'faces', d.surfFaces,'FaceVertexCData', d.baseRGB,'facecolor','interp','edgecolor','none','parent',d.hMainAxis);
    end
    
  case {d.hInnerAlphaEdit, d.hInnerAlphaSlider}
      set(d.hInnerSurface,'FaceAlpha',thisAlpha);
    
  case {d.hOuterAlphaEdit, d.hOuterAlphaSlider}
      set(d.hOuterSurface,'FaceAlpha',thisAlpha);
   
end
       
set(hFigure,'userdata',d);

%---------------------------------------------------- recomputes overlay object
function recomputeOverlay(handle,eventdata)

hFigure= get(handle,'parent');
d = get(hFigure,'userdata');

if any(d.alphaOverlayList)
  useAlpha = get(d.hOverlayAlphaCheck,'Value');
else
  useAlpha=false;
end
whichOverlays = get(d.hOverlayList,'Value');
overlayRoiModes = get(d.hOverlayRoiMode,'String');
overlayRoiMode = overlayRoiModes{get(d.hOverlayRoiMode,'Value')};

if any(d.alphaOverlayList) && (length(whichOverlays)>1 || ~d.alphaOverlayList(whichOverlays))
  set([d.hAlphaOverlayEditMin d.hAlphaOverlaySliderMin d.hAlphaOverlayEditMax d.hAlphaOverlaySliderMax],'enable','off');
end
if length(whichOverlays)>1
  set([d.hOverlayEditMin d.hOverlaySliderMin d.hOverlayEditMax d.hOverlaySliderMax],'enable','off')
else
% if only one overlay selected 
  %if not called from the one of the clipping controls,  allow clipping
  if ~ismember(handle,[d.hOverlayEditMin d.hOverlaySliderMin d.hOverlayEditMax d.hOverlaySliderMax])
    set([d.hOverlayEditMin d.hOverlaySliderMin d.hOverlayEditMax d.hOverlaySliderMax],'enable','on');
    set(d.hOverlayEditMin,'string',num2str(d.overlayCurrentClip(whichOverlays,1)),'enable','on');
    set(d.hOverlaySliderMin,'value',d.overlayCurrentClip(whichOverlays,1),'min',d.overlayClip(whichOverlays,1),'max',d.overlayClip(whichOverlays,2),'enable','on');
    set(d.hOverlayEditMax,'string',num2str(d.overlayCurrentClip(whichOverlays,2)),'enable','on');
    set(d.hOverlaySliderMax,'value',d.overlayCurrentClip(whichOverlays,2),'min',d.overlayClip(whichOverlays,1),'max',d.overlayClip(whichOverlays,2),'enable','on');
  else %otherwise, read value from the controls
    switch(handle)
      case d.hOverlaySliderMin
        d.overlayCurrentClip(whichOverlays,1) = get(handle,'value');
        set(d.hOverlayEditMin,'String',num2str(d.overlayCurrentClip(whichOverlays,1)));
      case d.hOverlaySliderMax
        d.overlayCurrentClip(whichOverlays,2) = get(handle,'value');
        set(d.hOverlayEditMax,'String',num2str(d.overlayCurrentClip(whichOverlays,2)));
      case d.hOverlayEditMin
        value = str2num(get(d.hOverlayEditMin,'string'));
        if  value<d.overlayClip(whichOverlays,1)
          set(d.hOverlayEditMin,'String',num2str(d.overlayClip(whichOverlays,1)));
          return;
        end
        d.overlayCurrentClip(whichOverlays,1) = value;
        set(d.hOverlaySliderMin,'Value',value);
      case d.hOverlayEditMax
        value = str2num(get(d.hOverlayEditMax,'string'));
        if  value>d.overlayClip(whichOverlays,2)
          set(d.hOverlayEditMin,'String',num2str(d.overlayClip(whichOverlays,2)));
          return;
        end
        d.overlayCurrentClip(whichOverlays,2) = value;
        set(d.hOverlaySliderMax,'value',value);
    end
    %and clip the data
    d.overlayClippedMask(:,whichOverlays) = d.overlayMaskData(:,whichOverlays);
    d.overlayClippedMask(d.overlayData(:,whichOverlays)<=d.overlayCurrentClip(whichOverlays,1)...
                       | d.overlayData(:,whichOverlays)>=d.overlayCurrentClip(whichOverlays,2),whichOverlays)=0;
  end
  if d.alphaOverlayList(whichOverlays) %same for alpha overlay clipping
    if ~ismember(handle,[d.hAlphaOverlayEditMin d.hAlphaOverlaySliderMin d.hAlphaOverlayEditMax d.hAlphaOverlaySliderMax])
    set([d.hAlphaOverlayEditMin d.hAlphaOverlaySliderMin d.hAlphaOverlayEditMax d.hAlphaOverlaySliderMax],'enable','on');
    set(d.hAlphaOverlayEditMin,'string',num2str(d.alphaOverlayCurrentClip(whichOverlays,1)),'enable','on');
    set(d.hAlphaOverlaySliderMin,'value',d.alphaOverlayCurrentClip(whichOverlays,1),'min',d.alphaOverlayClip(whichOverlays,1),'max',d.alphaOverlayClip(whichOverlays,2),'enable','on');
    set(d.hAlphaOverlayEditMax,'string',num2str(d.alphaOverlayCurrentClip(whichOverlays,2)),'enable','on');
    set(d.hAlphaOverlaySliderMax,'value',d.alphaOverlayCurrentClip(whichOverlays,2),'min',d.alphaOverlayClip(whichOverlays,1),'max',d.alphaOverlayClip(whichOverlays,2),'enable','on');
    else 
      set([d.hAlphaOverlayEditMin d.hAlphaOverlaySliderMin d.hAlphaOverlayEditMax d.hAlphaOverlaySliderMax],'enable','on');
      switch(handle)
        case d.hAlphaOverlaySliderMin
          d.alphaOverlayCurrentClip(whichOverlays,1) = get(handle,'value');
          set(d.hAlphaOverlayEditMin,'String',num2str(d.alphaOverlayCurrentClip(whichOverlays,1)));
        case d.hAlphaOverlaySliderMax
          d.alphaOverlayCurrentClip(whichOverlays,2) = get(handle,'value');
          set(d.hAlphaOverlayEditMax,'String',num2str(d.alphaOverlayCurrentClip(whichOverlays,2)));
        case d.hAlphaOverlayEditMin 
          value = str2num(get(d.hAlphaOverlayEditMin,'string'));
          if value<d.alphaOverlayClip(whichOverlays,1)
            set(d.hAlphaOverlayEditMin,'String',num2str(d.overlayClip(whichOverlays,1)));
            return;
          end
          d.alphaOverlayCurrentClip(whichOverlays,1) = value;
          set(d.hAlphaOverlaySliderMin,'Value',value);
        case d.hAlphaOverlayEditMax
          value = str2num(get(d.hAlphaOverlayEditMax,'string'));
          if value>d.alphaOverlayClip(whichOverlays,2)
            set(d.hAlphaOverlayEditMin,'String',num2str(d.overlayClip(whichOverlays,2)));
            return;
          end
          d.alphaOverlayCurrentClip(whichOverlays,2) = value;
          set(d.hAlphaOverlaySliderMax,'value',value);
      end
      d.overlayAlphaClippedMask(:,whichOverlays) = d.overlayAlphaMaskData(:,whichOverlays);
      d.overlayAlphaClippedMask(d.overlayAlphaData(:,whichOverlays)<=d.alphaOverlayCurrentClip(whichOverlays,1)...
                              | d.overlayAlphaData(:,whichOverlays)>=d.alphaOverlayCurrentClip(whichOverlays,2),whichOverlays)=0;
    end
  end
end

d = blendOverlays(d,whichOverlays,overlayRoiMode);

overlayAlpha = str2num(get(d.hOverlayAlphaEdit,'string'));
switch(handle)
  case d.hOverlayAlphaEdit
    overlayAlpha = str2num(get(handle,'string'));
     set(d.hOverlayAlphaSlider,'value',overlayAlpha);
  case d.hOverlayAlphaSlider
    overlayAlpha = get(handle,'value');
    set(d.hOverlayAlphaEdit,'string',num2str(overlayAlpha));
end

if useAlpha
  cubesAlphaData = overlayAlpha*d.cubesAlphaData;
else
  cubesAlphaData = overlayAlpha*ones(size(d.cubesAlphaData));
end
    
if ~isempty(cubesAlphaData)
   [faceXcoords, faceYcoords, faceZcoords,faceCData,faceAlphaData] = ...
         makeCubeFaces(d.cubesScanCoords,d.cubesRGB',cubesAlphaData',all(cubesAlphaData(:)==1)||d.baseType==1,d.base2scan);
   %to use transparency, I need to convert X,Y,Z data to face/vertex data
   %patchStruct = surf2patch(faceXcoords,faceYcoords,faceZcoords,faceCData); %surf2patch does not handle color data correctly
   patchStruct.vertices = [reshape(faceXcoords,numel(faceXcoords),1) reshape(faceYcoords,numel(faceYcoords),1) reshape(faceZcoords,numel(faceZcoords),1)];
   patchStruct.faces = reshape((1:size(patchStruct.vertices,1)),[4 size(faceXcoords,2)])';
else
   patchStruct.vertices = [];
   patchStruct.faces = [];
   faceCData = 1;
   faceAlphaData = 1;
end

if isfield(d,'hOverlay')
   set(d.hOverlay,'vertices',patchStruct.vertices,'faces',patchStruct.faces,...
     'FaceVertexCData',faceCData','FaceVertexAlphaData',faceAlphaData','FaceAlpha','flat','faceColor','flat');
else
   patchStruct.facevertexcdata = faceCData';
   d.hOverlay = patch(patchStruct,'edgecolor','none','FaceVertexAlphaData',faceAlphaData','FaceAlpha','flat','faceColor','flat','parent',d.hMainAxis);
end
      
% %       %count number of visible voxels in each ROI
% %       for iRoi = 1:length(d.roiSize)
% %          set(d.hRoiVoxelNum(iRoi),'String',[num2str(sum(d.overlayRoiIndex{iRoi}>=firstIndex & d.overlayRoiIndex{iRoi}<=lastIndex)) '/' num2str(d.roiSize(iRoi))]);
% %       end

set(hFigure,'userdata',d);

return;

%--------------------------------------returns the coordinates of the faces of cubes centered on voxel coordinates
function d = blendOverlays(d,whichOverlays,overlayRoiMode)

switch(overlayRoiMode)
  case 'Show all voxels'
    overlayAlphaMaskData = d.overlayAlphaClippedMask;
    overlayAlpha = d.overlayAlpha;
    overlayMaskData = d.overlayClippedMask;
%     overlayData = d.overlayData;
    overlayRGB = d.overlayRGB;
    overlayScanCoords = d.overlayScanCoords;
    
  case 'Restrict to all ROIs'
    overlayAlphaMaskData = d.overlayAlphaClippedMask(d.overlayAllRoisIndex,:);
    overlayAlpha = d.overlayAlpha(d.overlayAllRoisIndex,:);
    overlayMaskData = d.overlayClippedMask(d.overlayAllRoisIndex,:);
%     overlayData = d.overlayData(d.overlayAllRoisIndex,:);
    overlayRGB = d.overlayRGB(d.overlayAllRoisIndex,:,:);
    overlayScanCoords = d.overlayScanCoords(:,d.overlayAllRoisIndex);
    
  case 'Restrict to selected ROIs'
    overlayRoisIndex = [];
    for iRoi=get(d.hRoiList,'value')
      overlayRoisIndex = union(d.overlayRoiIndex{iRoi},overlayRoisIndex);
    end
    overlayAlphaMaskData = d.overlayAlphaClippedMask(overlayRoisIndex,:);
    overlayAlpha = d.overlayAlpha(overlayRoisIndex,:);
    overlayMaskData = d.overlayClippedMask(overlayRoisIndex,:);
%     overlayData = d.overlayData(overlayRoisIndex,:);
    overlayRGB = d.overlayRGB(overlayRoisIndex,:,:);
    overlayScanCoords = d.overlayScanCoords(:,overlayRoisIndex);
    
end

%%%%%%%%%%%%%%%%%%% Compute cube faces

%apply alpha mask to alpha overlay
overlayAlpha(~overlayAlphaMaskData)=0;

cOverlay = 0;
for iOverlay=whichOverlays
  cOverlay = cOverlay+1;
    thisAlpha = overlayAlpha(:,iOverlay);
    % apply the mask to alpha data
    thisAlpha(~overlayMaskData(:,iOverlay))=0;
  % 1) pre-multiply colors by alpha
    RGB = repmat(thisAlpha,[1 3]).*overlayRGB(:,:,iOverlay);
    if cOverlay==1
      d.cubesRGB=RGB;
      d.cubesAlphaData = thisAlpha;
    else
    % 2) add pre-multiplied colormaps (and alpha channels)
      d.cubesRGB = d.cubesRGB.*(1-RGB)+ RGB;
      d.cubesAlphaData = d.cubesAlphaData .* (1-thisAlpha)+thisAlpha;
    end
end

%here we divide by the computed alpha, because there is no base to multiply it with
d.cubesRGB = d.cubesRGB./repmat(d.cubesAlphaData,[1 3]);
%exclude Nan overlay values
voxelsToKeepIndices = find(any(~isnan(d.cubesRGB),2));
%sort according to alpha values, this will be important when removing duplicate cube faces
d.cubesAlphaData = d.cubesAlphaData(voxelsToKeepIndices);
[d.cubesAlphaData, sortIndex] = sort(d.cubesAlphaData);
voxelsToKeepIndices = voxelsToKeepIndices(sortIndex);

%apply excluding and/or sorting to other RGB and alpha data arrays
d.cubesRGB = d.cubesRGB(voxelsToKeepIndices,:);
d.cubesScanCoords = overlayScanCoords(:,voxelsToKeepIndices);

%--------------------------------------returns the coordinates of the faces of cubes centered on voxel coordinates
function [facesXcoords, facesYcoords, facesZcoords, colorData, alphaData] = makeCubeFaces(voxelCoords,colorData,alphaData,opaque,base2scan)

%coordinates of a cube around each remaining voxel
cubeX = [-.5 -.5 .5 .5;-.5 -.5 .5 .5;-.5 -.5 -.5 -.5;.5 .5 .5 .5;-.5 .5 .5 -.5;-.5 .5 .5 -.5]';
cubeY = [-.5 -.5 -.5 -.5;.5 .5 .5 .5;-.5 .5 .5 -.5;-.5 .5 .5 -.5;-.5 -.5 .5 .5;-.5 -.5 .5 .5]';
cubeZ = [-.5 .5 .5 -.5;-.5 .5 .5 -.5;-.5 -.5 .5 .5;-.5 -.5 .5 .5;-.5 -.5 -.5 -.5;.5 .5 .5 .5]';

voxelsNumber = size(voxelCoords,2);
cubesXcoords = repmat(shiftdim(voxelCoords(1,:),-1),[4 6 1]) + repmat(cubeX,[1 1 voxelsNumber]);
cubesYcoords = repmat(shiftdim(voxelCoords(2,:),-1),[4 6 1]) + repmat(cubeY,[1 1 voxelsNumber]);
cubesZcoords = repmat(shiftdim(voxelCoords(3,:),-1),[4 6 1]) + repmat(cubeZ,[1 1 voxelsNumber]);

facesXcoords = reshape(cubesXcoords,4,6*voxelsNumber);
facesYcoords = reshape(cubesYcoords,4,6*voxelsNumber);
facesZcoords = reshape(cubesZcoords,4,6*voxelsNumber);

if ~ieNotDefined('colorData')
  colorData = reshape(repmat(colorData,[6 1]),size(colorData,1),6*voxelsNumber);
end
if ~ieNotDefined('alphaData')
  alphaData = reshape(repmat(alphaData,[6 1]),1,6*voxelsNumber);
end

%find duplicate faces to remove (one of each duplicate)
facesCoords = [facesXcoords; facesYcoords; facesZcoords]';
[uniqueFaceCoords, uniqueIndices] = unique(facesCoords,'rows');

if opaque   %if everything is opaque, then we can just keep the outside faces, that is remove the other set of duplicate faces
   doubleIndices = setdiff(1:size(facesCoords,1),uniqueIndices);
   [uniqueFaceCoords, uniqueIndices] = setdiff(facesCoords,facesCoords(doubleIndices,:),'rows');
end
if ~ieNotDefined('colorData')
  colorData = colorData(:,uniqueIndices);
end
if ~ieNotDefined('alphaData')
  alphaData = alphaData(uniqueIndices);
  %round any alpha value >.985 to 1, because of a bug in the display
  alphaData(alphaData>.985)=1;
else
  alphaData = ones(size(uniqueIndices))';
end

if ~isempty(uniqueFaceCoords)
  facesXcoords = uniqueFaceCoords(:,1:4)';
  facesYcoords = uniqueFaceCoords(:,5:8)';
  facesZcoords = uniqueFaceCoords(:,9:12)';
end

[facesXcoords,facesYcoords,facesZcoords] = convert2baseCoords(facesXcoords,facesYcoords,facesZcoords,base2scan);

%-----------------------
function [baseX,baseY,baseZ] = convert2baseCoords(X,Y,Z,base2scan)

baseX = reshape(X,1,numel(X));
baseY = reshape(Y,1,numel(Y));
baseZ = reshape(Z,1,numel(Z));

baseCoords = base2scan\[baseX;baseY;baseZ;ones(1,size(baseX,2))];

baseX = baseCoords(1,:);
baseY = baseCoords(2,:);
baseZ = baseCoords(3,:);

baseX = reshape(baseX,size(X));
baseY = reshape(baseY,size(Y));
baseZ = reshape(baseZ,size(Z));


