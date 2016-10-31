% correlationMap.m
%
%      usage: makeCorrelationMap()
%         by: justin gardner
%       date: 05/12/08
%    purpose: interrogator function that creates a correlation map from
%             the voxel clicked on (correlation is based on event related
%             responses).
%
function retval = makeCorrelationMap(v,overlayNum,scan,x,y,s,roi,varargin)

% get the view number
viewNum = viewGet(v,'viewNum');

% if we have it
if ismember(viewGet(v,'analysisType'),{'erAnal','deconvAnal');
  
  % then get the d strucutre
  analysis = viewGet(v,'analysis',viewGet(v,'curAnalysis'));
  d = analysis.d{scan};
  if isempty(d)
    disp(sprintf('(eventRelatedPlot) Event related not for scan %i',scan));
    return
  end
  d.r2 = analysis.overlays(1).data{scan};
  d.dim = viewGet(v,'scanDims');
  roiName = '';roiN = 0;
  
  % now, pull out hdr for the current voxel (note that this is
  % an interlaced trace -- i.e. the first time point is the first hdr
  % the second time point is the second hdr etc.) This makes it
  % so that it will be the same as the targetHDR reshaped below.
  % Also, it won't matter for correlation what the order in the array is
  % as long as it matches the targetHDR
  if isempty(roi)
    sourceHDR = reshape(squeeze(d.ehdr(x,y,s,:,:)),1,d.nhdr*d.hdrlen);
  else
    % keep name of roi
    roiName = roi{1}.name;
    % for this case, we get all the ehdrs for the roi
    disppercent(-inf,sprintf('(makeCorrelationMap) Using roi %s as source',roiName));
    roi{1}.scanCoords = getROICoordinates(v,roi{1});
    roiN = size(roi{1}.scanCoords,2); 
    for i = 1:roiN
      xi = roi{1}.scanCoords(1,i);
      yi = roi{1}.scanCoords(2,i);
      si = roi{1}.scanCoords(3,i);
      sourceHDR(i,:) = reshape(squeeze(d.ehdr(xi,yi,si,:,:)),1,d.nhdr*d.hdrlen);
      disppercent(i/roiN);
    end
    disppercent(inf);
  end
  
  % now reshape the ehdr matrix as a matrix where the first dimension
  % is the voxel number and the second dimension is the response
  targetHDR = reshape(d.ehdr,prod(d.dim(1:3)),d.nhdr*d.hdrlen);

  % compute correlation using projectOutMeanVector
  p = projectOutMeanVector([],sourceHDR,targetHDR);

  % reshape the ehdr's that have the projected out one
  d.projectedOutEhdr = reshape(p.tSeries,d.dim(1),d.dim(2),d.dim(3),d.nhdr,d.hdrlen);

  % get the projection ehdr
  d.projectionEhdr = reshape(p.sourceMeanVector,d.nhdr,d.hdrlen);

  % reshape the correlation map
  r = reshape(p.r,d.dim(1),d.dim(2),d.dim(3));
  
  % install the overlay
  mrDispOverlay(r,viewGet(v,'curScan'),viewGet(v,'curGroup'),v,'overlayName=r','interrogator=erCorrelationPlot','params.x',x,'params.y',y,'params.s',s,'params.roiName',roiName,'params.roiN',roiN,'range',[-1 1],'clip',[0.1 -0.1],'cmap',[flipud(fliplr(hot(128)));hot(128)],'colormapType=setRangeToMaxAroundZero','d',d);

  % refresh display
  refreshMLRDisplay(viewGet(v,'viewNum'));
 
end

