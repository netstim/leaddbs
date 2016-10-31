% computeOverlay.m
%
%        $Id$
%      usage: overlays = computeOverlay(thisView,base2overlay,baseCoordsHomogeneous,baseDims,<alpha>)
%         by: David Heeger
%       date: 10/16/07
%    purpose: this used to live within the refreshMLRDisplay
%             function, but has been pulled out so that it
%             can be called by other functions (e.g. mrPrint)-jg
%
function overlays = computeOverlay(thisView,base2overlay,baseCoordsHomogeneous,baseDims,alpha)

% check arguments
if ~any(nargin == [4 5])
  help computeOverlay
  return
end

% viewGet some stuff
curOverlays = viewGet(thisView,'currentOverlay');
nCurOverlays = length(curOverlays);
analysisNum = viewGet(thisView,'currentAnalysis');
scan = viewGet(thisView,'curscan');

sliceInfo.base2overlay = base2overlay;
sliceInfo.baseCoordsHomogeneous = baseCoordsHomogeneous;
sliceInfo.baseDims = [baseDims 1];

if ~isempty(curOverlays) && scan && viewGet(thisView,'nscans') % there needs to be at least one scan in the group
  curAlphaOverlays = zeros(1,nCurOverlays);
  for iOverlay = 1:nCurOverlays
    thisNum = viewGet(thisView,'overlayNum',viewGet(thisView,'alphaOverlay',curOverlays(iOverlay)));
    if ~isempty(thisNum)
      curAlphaOverlays(iOverlay) = thisNum;
    end
  end
  curOverlays = [curOverlays curAlphaOverlays];
  
  if strcmp(viewGet(thisView,'overlayInterpFunction',analysisNum),'corAnalInterp') ...    %if it's a correlation analysis
    && size(sliceInfo.baseCoordsHomogeneous,3)>1 ...                                    %and there are several depths
    && strcmp(mrGetPref('multiSliceProjectionMethod'),'Average')                        %and we need to average across depths
    %check if overlays are phase or amplitude    
    corAnalOverlays=zeros(size(curOverlays));
    cOverlay=0;
    for iOverlay=curOverlays
      cOverlay = cOverlay+1;
      if iOverlay
        overlayType{cOverlay} = viewGet(thisView,'overlayType',iOverlay,analysisNum);
        % we assume that phase overlays always follow amplitude overlays
        switch(overlayType{cOverlay})
          case 'amp'
            corAnalOverlays(cOverlay) = curOverlays(cOverlay)+1;
          case 'ph'
            corAnalOverlays(cOverlay) = curOverlays(cOverlay)-1;
        end
      else
        overlayType{cOverlay}='none';
      end
    end
    curOverlays = [curOverlays corAnalOverlays];
  end

  %get overlays and masks of clipped values 
  [overlayMasks,overlayImages,overlays.coords]= maskOverlay(thisView,curOverlays,scan,sliceInfo);
else
  overlayImages = [];
end

if ~isempty(overlayImages)
  overlays.coords = overlays.coords{1};
  alphaOverlayMasks = overlayMasks{1}(:,:,:,nCurOverlays+1:2*nCurOverlays);
  overlayMasks = overlayMasks{1}(:,:,:,1:nCurOverlays);
  
  if size(overlayImages{1},3)>1
    overlayMasks = any(overlayMasks,3);
    alphaOverlayMasks = any(alphaOverlayMasks,3);
    switch(mrGetPref('multiSliceProjectionMethod'))
      case 'Average'
        if strcmp(viewGet(thisView,'overlayInterpFunction',analysisNum),'corAnalInterp')
          corAnalOverlayImages = overlayImages{1}(:,:,:,2*nCurOverlays+1:4*nCurOverlays);
          overlayImages = overlayImages{1}(:,:,:,1:2*nCurOverlays);
          %for amplitude and phase in a correlation analysis, transform into complex numbers
          for iOverlay=1:2*nCurOverlays
            if corAnalOverlays(iOverlay) %if we have the other part of the complex data (phase or amplitude)
              switch(overlayType{iOverlay})
                case 'amp'
                  overlayImages(:,:,:,iOverlay) = overlayImages(:,:,:,iOverlay) .* exp(1i*corAnalOverlayImages(:,:,:,iOverlay));
                case 'ph'
                  overlayImages(:,:,:,iOverlay) = corAnalOverlayImages(:,:,:,iOverlay) .* exp(1i*overlayImages(:,:,:,iOverlay));
              end
            elseif strcmp(overlayType{iOverlay},'ph') %if we don't have the amplitude data but it's a phase overlay
              overlayImages(:,:,:,iOverlay) = exp(1i*overlayImages(:,:,:,iOverlay));
            end
          end
        else
          overlayImages = overlayImages{1};
        end
        overlayImages = cast(nanmean(double(overlayImages),3),class(overlayImages));
        if strcmp(viewGet(thisView,'overlayInterpFunction',analysisNum),'corAnalInterp')
          for iOverlay=1:nCurOverlays
            switch(overlayType{iOverlay})
              case 'amp'
                overlayImages(:,:,:,iOverlay) = abs(overlayImages(:,:,:,iOverlay)); %keep only the amplitude
              case 'ph'
                ang = angle(overlayImages(:,:,:,iOverlay));%keep only the phase
                ang(ang < 0) = ang(ang < 0)+2*pi;
                overlayImages(:,:,:,iOverlay) = ang;
            end
          end
        end
        alphaOverlayImages = overlayImages(:,:,:,nCurOverlays+1:2*nCurOverlays);
        overlayImages = overlayImages(:,:,:,1:nCurOverlays);
        
      case 'Maximum Intensity Projection'
        %maybe could have a special case for phase data (in corAnal)
        alphaOverlayImages = overlayImages{1}(:,:,:,nCurOverlays+1:2*nCurOverlays);
        overlayImages = overlayImages{1}(:,:,:,1:nCurOverlays);
        alphaOverlayImages= permute(alphaOverlayImages,[4 3 1 2]);
        newAlphaOverlayImages = NaN([size(overlayImages,4) size(overlayImages,1) size(overlayImages,2)]);
        newOverlayImages = NaN([size(overlayImages,1) size(overlayImages,2) 1 size(overlayImages,4)]);
        for iOverlay = 1:nCurOverlays
          [newOverlayImages(:,:,1,iOverlay),maxIndex] = max(overlayImages(:,:,:,iOverlay),[],3);
          %get the corresponding alpha values (?)
          for iSlice = 1:size(overlayImages,3)
            thisIndex = find(maxIndex==iSlice);
            newAlphaOverlayImages(iOverlay,thisIndex) = alphaOverlayImages(iOverlay,iSlice,thisIndex);
          end
        end
        alphaOverlayImages = permute(newAlphaOverlayImages,[2 3 4 1]);
        overlayImages = newOverlayImages;
    end
  else
    alphaOverlayImages = overlayImages{1}(:,:,:,nCurOverlays+1:2*nCurOverlays);
    overlayImages = overlayImages{1}(:,:,:,1:nCurOverlays);
  end
  
  overlays.overlayIm = permute(overlayImages,[1 2 4 3]); 
  imageDims = [size(overlayImages,1) size(overlayImages,2)];
  
  
  %alpha maps based on each overlay's alphaoverlay or alpha value
  overlays.alphaMaps=repmat(NaN([imageDims 1 nCurOverlays],mrGetPref('defaultPrecision')),[1 1 3 1]);
  overlays.RGB=repmat(NaN([imageDims 1 nCurOverlays],mrGetPref('defaultPrecision')),[1 1 3 1]);
  
  % Finally, make the alphaMap. Normally this is just set
  % to 0 or 1 so that voxels are hard thresholded. If the
  % overlays has an alphaOverlay field set to the name
  % of another overlays, then we will use the values in
  % that overlays to set the alpha
  alphaOverlayExponent = viewGet(thisView,'alphaOverlayExponent');
  for iOverlay = 1:nCurOverlays
    if nargin == 4
      alpha = viewGet(thisView,'alpha',curOverlays(iOverlay));
    end
    if ~curAlphaOverlays(iOverlay)
      overlays.alphaMaps(:,:,:,iOverlay) = repmat(alpha*overlayMasks(:,:,:,iOverlay),[1 1 3]);
    else
      %get the alpha overlay data
      alphaOverlayImage = alphaOverlayImages(:,:,:,iOverlay);
      % get the range of the alpha overlays
      range = viewGet(thisView,'overlayRange',curAlphaOverlays(iOverlay));
      % handle setRangeToMax
      if strcmp(viewGet(thisView,'overlayCtype',curAlphaOverlays(iOverlay)),'setRangeToMax')
        maxRange = max(clip(1),min(alphaOverlayImage(overlayMasks(:,:,:,iOverlay))));
        if ~isempty(maxRange),range(1) = maxRange;end
        minRange = min(max(alphaOverlayImage(overlayMasks(:,:,:,iOverlay))),clip(2));
        if ~isempty(minRange),range(2) = minRange;end
      end
      % now compute the alphaOverlay as a number from
      % 0 to 1 of the range
      alphaOverlayImage = ((alphaOverlayImage-range(1))./diff(range));
      if alphaOverlayExponent(iOverlay)<0   % if the alpha overlays exponent is negative, set it positive and take the 1-alpha for the alpha map
         alphaOverlayExponent(iOverlay) = -alphaOverlayExponent(iOverlay);
         alphaOverlayImage = 1-alphaOverlayImage;
      end
      alphaOverlayImage = alpha*(alphaOverlayImage.^alphaOverlayExponent(iOverlay));
      alphaOverlayImage(isnan(alphaOverlayImage)) = 0;
      alphaOverlayImage(alphaOverlayImage>alpha) = alpha;
      alphaOverlayImage(alphaOverlayImage<0) = 0;
      %mask the alpha overlays with its own mask and the overlay mask
      alphaOverlayImage = alphaOverlayImage.*overlayMasks(:,:,:,iOverlay);
      overlays.alphaMaps(:,:,:,iOverlay) = repmat(alphaOverlayImage.*alphaOverlayMasks(:,:,:,iOverlay),[1 1 3]); 
    end   

    % Rescale current overlay.
    overlays.cmap(:,:,iOverlay) = viewGet(thisView,'overlayCmap',curOverlays(iOverlay));
    overlayCType = viewGet(thisView,'overlayCtype',curOverlays(iOverlay));
    % For multiAxis we need to have a consistent colorbar, so we force
    % setRangeToMax to go across slices
    if (viewGet(thisView,'baseType') == 0) && viewGet(thisView,'baseMultiAxis') && isequal(overlayCType,'setRangeToMax')
      overlayCType = 'setRangeToMaxAcrossSlices';
    end
    switch(overlayCType)
      case 'normal'
        overlays.range(iOverlay,:) = viewGet(thisView,'overlayRange',curOverlays(iOverlay));
      case 'setRangeToMax'
        overlays.range(iOverlay,:) = viewGet(thisView,'overlayClip',curOverlays(iOverlay));
        if ~isempty(overlays.overlayIm(overlayMasks(:,:,:,iOverlay)))
          thisOverlayIm=overlays.overlayIm(:,:,iOverlay);
          overlays.range(iOverlay,1) = max(overlays.range(iOverlay,1),min(thisOverlayIm(overlayMasks(:,:,:,iOverlay))));
          overlays.range(iOverlay,2) = min(max(thisOverlayIm(overlayMasks(:,:,:,iOverlay))),overlays.range(iOverlay,2));
        end
      case 'setRangeToMaxAroundZero'
       if (viewGet(thisView,'baseType') == 0) && viewGet(thisView,'baseMultiAxis')
	 oneTimeWarning('setRangeToMaxWithMultiAxis','(computeOverlay) !!! Overlay set to setRangeToMaxAroundZero and you are displaying with multiple axis. Note that the colorbar will only be accurate for the axial view. Consider setting Edit/Overlay/Edit OVerlay to setRangeToMaxAcrossSlices. !!!',1);
       end
        overlays.range(iOverlay,:) = viewGet(thisView,'overlayClip',curOverlays(iOverlay));
        if ~isempty(overlays.overlayIm(overlayMasks(:,:,:,iOverlay)))
          maxval = max(abs(overlays.overlayIm(overlayMasks(:,:,:,iOverlay))));
          overlays.range(iOverlay,1) = -maxval;
          overlays.range(iOverlay,2) = maxval;
        end
      case 'setRangeToMaxAcrossSlices' %this doesn't take into account voxels that would be masked in other slices than the one displayed
        overlays.range(iOverlay,:) = viewGet(thisView,'overlayClip',curOverlays(iOverlay));
        minOverlay = viewGet(thisView,'minOverlaydata',curOverlays(iOverlay),analysisNum,scan); 
        maxOverlay = viewGet(thisView,'maxOverlaydata',curOverlays(iOverlay),analysisNum,scan);
        if ~isempty(minOverlay)
          overlays.range(iOverlay,1) = max(minOverlay,overlays.range(1));
        end
        if ~isempty(maxOverlay)
          overlays.range(iOverlay,2) = min(maxOverlay,overlays.range(2));
        end
      case 'setRangeToMaxAcrossSlicesAndScans' %same here
        overlays.range(iOverlay,:) = viewGet(thisView,'overlayClip',curOverlays(iOverlay));
        minOverlay = viewGet(thisView,'minOverlaydata',curOverlays(iOverlay),analysisNum); 
        maxOverlay = viewGet(thisView,'maxOverlaydata',curOverlays(iOverlay),analysisNum);
        if ~isempty(minOverlay)
          overlays.range(iOverlay,1) = max(minOverlay,overlays.range(1));
        end
        if ~isempty(maxOverlay)
          overlays.range(iOverlay,2) = min(maxOverlay,overlays.range(2));
        end
    end
    overlays.RGB(:,:,:,iOverlay) = rescale2rgb(overlayImages(:,:,:,iOverlay),overlays.cmap(:,:,iOverlay),overlays.range(iOverlay,:));
      
  end
  %alpha based on the alpha slider
  if nargin == 4
    alpha = viewGet(thisView,'alpha');
  end
  overlays.alphaMap = alpha*any(overlays.alphaMaps,4);
  %overlays.alphaMap = repmat(alpha*ones(imageDims,mrGetPref('defaultPrecision')),[1 1 3]); 

else
  overlays.coords = [];
  overlays.overlayIm = [];
  overlays.RGB = [];
end

