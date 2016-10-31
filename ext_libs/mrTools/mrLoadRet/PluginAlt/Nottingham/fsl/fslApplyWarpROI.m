% fslApplyWarpROI.m
%
%        $Id: fslApplyWarpROI.m 2119 2011-05-09 22:46:58Z julien $ 
%      usage: fslApplyWarpROI(thisView)
%         by: julien besle
%       date: 06/10/2010
%    purpose: applies non-linear FSL registration to chosen ROI(s)
%
function thisView = fslApplyWarpROI(thisView)

if strcmp(mrGetPref('fslPath'),'FSL not installed')
  mrWarnDlg('(fslApplyWarpROI) No path was provided for FSL. Please set MR preference ''fslPath'' by running mrSetPref(''fslPath'',''yourpath'')')
  return
end

useApplyWarp = 1; %The solution I chose is to create a volume of zeros and put ones where the ROI is, 
%then run this through applywarp using the spline coefficient  and nearest neighbour interpolation
% and finally find the coordinates of the transformed ROI (where the ones are)
% if useApplyWarp=0, coordinates are converted directly using warp fields created using fnirtfileutils
% but this does not preserve the volume of the ROI and takes more time.

currentROIName = viewGet(thisView,'roiname');
needToRefresh=0;

numBases = viewGet(thisView,'numberofBaseVolumes');
for iBase = 1:numBases
   baseNames{iBase} = viewGet(thisView,'baseName',iBase);
end
baseNames = [{'current scan'},baseNames];

keepAsking = 1;
while keepAsking
  
  params = {...
   {'fnirtInputSpace',baseNames,'type=popupmenu','A volume in a space equivalent to the volume that served as the FNIRT input volume, generally the scan space.'},...
   {'destinationSpace',{'FNIRT input space','Keep original ROI space'},'type=popupmenu','Which coordinate space the transformed ROI will be in. If the destination space is different from the fnirt input space, the total ROI volume will be conserved, but results might not be as accurate.'},...
   {'roiNamePrefix','FNIRT_','String to append at the start of each saved roi name'},...
   {'roiNameSuffix','','String to append at the end of each saved roi name'},...
          };
  params = mrParamsDialog(params, 'fslApplyWarpROI parameters');
  % Abort if params empty
  if ieNotDefined('params')
    return
  else
  
%   baseNum = [1 2];
%   while length(baseNum)>1
%     baseNum = selectInList(thisView,'bases','Select the input volume of the FNIRT warp coeffs (or equivalent)');
%     if length(baseNum)>1
%       mrWarnDlg('(fslApplyWarpROI) Please select only one base');
%     end
%   end
%   if isempty(baseNum)
%     return;
%   else

    while keepAsking
      roiList = selectInList(thisView,'rois','Select ROI(s) to warp');
      if isempty(roiList)
         break;
      else
         keepAsking = 0;
         [warpCoefFilename warpCoefPathname] = uigetfile('*.img;*.nii','FNIRT Warp spline coefficient file');
         if isnumeric(warpCoefFilename)
            keepAsking=1;
         end
      end
    end
  end
end

%find temporary file extension based on FSL preference
switch(getenv('FSLOUTPUTTYPE'))
  case 'NIFTI'
    tempFilename='temp.nii';
  case 'NIFTI_PAIR'
    tempFilename='temp.img';
  case ''
    mrWarnDlg('(fslApplyWarpSurfOFF) Environment variable FSLOUTPUTTYPE is not set');
    return;
  otherwise
    mrWarnDlg(sprintf('(fslApplyWarpOverlays) Environment variable FSLOUTPUTTYPE is set to an unsupported value (%s)',getenv('FSLOUTPUTTYPE')));
    return;
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%% using warping coefficient and applywarp
[~,baseNum]=ismember(params.destinationSpace,baseNames);
baseNum = baseNum-1;
if baseNum>0
  fnirtInputXform = viewGet(thisView,'basexform',baseNum);
  fnirtInputDims = viewGet(thisView,'basedims',baseNum); 
  fnirtInputHdr = viewGet(thisView,'basehdr',baseNum);
  fnirtInputVoxelSize = viewGet(thisView,'baseVoxelSize',baseNum);
else
  scanNum = viewGet(thisView,'curscan');
  fnirtInputXform = viewGet(thisView,'scanxform',scanNum);
  fnirtInputDims = viewGet(thisView,'scandims',scanNum); 
  fnirtInputHdr = viewGet(thisView,'niftihdr',scanNum);
  fnirtInputVoxelSize = viewGet(thisView,'scanVoxelSize',scanNum);
end

fprintf('\n');
for iRoi = 1:length(roiList)
  thisRoi = viewGet(thisView,'ROI',roiList(iRoi));
    
  %transformation matrix between roi space and fnirt input volume space
  roi2fnirtInput = fnirtInputXform\thisRoi.xform;
    
  %for some reason, there might be NaNs in the roi coords, remove those
  thisRoi.coords = thisRoi.coords(:,~any(isnan(thisRoi.coords(1:3,:)),1));
  if isempty(thisRoi.coords)
    mrwarnDlg(['Roi ' thisRoi.name ' is empty, skipping ...']);
    fprintf('\n');
  else
    thisWarpedRoi = thisRoi;
    thisWarpedRoi.name = [params.roiNamePrefix thisWarpedRoi.name params.roiNameSuffix];
    baseCoords = round(roi2fnirtInput*thisRoi.coords);
    outsideVoxels = any(baseCoords(1:3,:)< ones(3,size(baseCoords,2)) |...
    baseCoords(1:3,:)> repmat(fnirtInputDims',1,size(baseCoords,2)));
    if any(outsideVoxels) 
      mrWarnDlg(['(fslApplyWarpROI) ' num2str(nnz(outsideVoxels)) ' voxels in ' thisRoi.name ' are outside the input FNIRT volume, removing...']);
      thisRoi.coords(:,outsideVoxels)=[];
    end
    if useApplyWarp
      if any(any( (roi2fnirtInput - eye(4))>1e-6)) && strcmp(params.destinationSpace,'Keep original ROI space')
        disp('(fslApplyWarpROI) ROI and FNIRT input volume spaces are different. Using partial volume correction...')
        %here we will be following the logic of xformROIcoords, supersampling the roi space, transforming
        % each sub-voxel volume into the FNIRT input volume and estimating the accumulated partial voluming at each point 
        % the difference will be that we will FNIRT the continous accumulated values, and then apply a threshold that conserve the total ROI volume
        superSampling = [5 5 5];
        centeredSamples = floor(superSampling/2)./superSampling;%assumes that supersampling factors are odd
        sampleSizes = 1./superSampling;
        roiVolume = zeros(fnirtInputDims,'single');
        % for each subsample of all roi coordinates, compute the coordinates in the FNIRT input volume 
        % and count how many times any input coordinate falls into any voxel of this volume (accumulated volume)
        hWait = mrWaitBar(-inf,['(fslApplyWarpROI) Computing ROI partial voluming (supersampling = ' mat2str(superSampling) ')']);
        counter=0;
        for i=-centeredSamples(1):sampleSizes(1):centeredSamples(1)
          for j=-centeredSamples(2):sampleSizes(2):centeredSamples(2)
            for k=-centeredSamples(3):sampleSizes(3):centeredSamples(3)
              counter=counter+1;
              baseCoords = round(roi2fnirtInput*(thisRoi.coords+repmat([i j k 0]',1,size(thisRoi.coords,2))));
              baseIndices = mrSub2ind(fnirtInputDims(1:3),baseCoords(1,:),baseCoords(2,:),baseCoords(3,:)); %use mrSub2ind to identify coordinates outside the scan
              baseIndices = baseIndices(~isnan(baseIndices)); 
              for jj=1:length(baseIndices)
                roiVolume(baseIndices(jj)) = roiVolume(baseIndices(jj)) + 1;
              end
              mrWaitBar(counter/prod(superSampling),hWait);
            end
          end
        end
        mrCloseDlg(hWait)
        %apply FNIRT warp to the continuous accumulated values
        warpedRoiVolume = fslApplyWarp(roiVolume, [warpCoefPathname warpCoefFilename], tempFilename, fnirtInputHdr, 'sinc');

        %now interpolate into the roi volume (not sure this would be accurate if voxels are not isometric..)
        % find ROI coordinates bounding box in roi space
        minRoiCoords = min(thisRoi.coords(1:3,:),[],2);
        maxRoiCoords = max(max(thisRoi.coords(1:3,:),[],2),minRoiCoords+1); %(make sure the box is at least 2*2*2)
        boxDims = (maxRoiCoords - minRoiCoords +1)';
        [X,Y,Z]=ndgrid(minRoiCoords(1):maxRoiCoords(1),minRoiCoords(2):maxRoiCoords(2),minRoiCoords(3):maxRoiCoords(3));
        %transform it to base space
        boxBaseCoords=roi2fnirtInput*[X(:)';Y(:)';Z(:)';ones(1,numel(X))];
        boxVolume = interp3(warpedRoiVolume,boxBaseCoords(2,:),boxBaseCoords(1,:),boxBaseCoords(3,:),'*cubic');
        boxVolume(isnan(boxVolume))=0;

        %threshold the ROI in order to conserve the total ROI volume (=number of voxels since voxel size didn't change),
        %removing voxels with the least volume accumulation
        %this could be a bit of a problem if the non-linear transformation kicks a lot of ROI voxels out of the input volume
        % so we'll multiply by a correcting factor
        roiSize = round(size(thisRoi.coords,2)*sum(warpedRoiVolume(:))/sum(roiVolume(:)));
        if roiSize>nnz(boxVolume) %I don't think that would happen, unless FNIRT moves a lot of voxels out of the input volume
          disp('(fslApplyWarpROI) Not enough partial volume voxels for new roi size, entering debug mode');
          keyboard;
        end
        [~,indices] = sort(boxVolume,'descend');
  %       arrayviewer(reshape(boxVolume,boxDims));
        [x,y,z] = ind2sub(boxDims,indices(1:roiSize));
        thisWarpedRoi.coords = [minRoiCoords(1)-1+x;minRoiCoords(2)-1+y;minRoiCoords(3)-1+z;ones(1,roiSize)];

  %       thresholdBox = false(boxDims);
  %       thresholdBox(indices(1:size(thisRoi.coords,2)))=true;
  %       arrayviewer(thresholdBox);
  
      else %if everything is in the same space or the roi is to be converted to the input FNIRT space anyway
        %we simply use applyWarp with nearest neighbor interpolation, which seems to give better results in most cases
        %we require that ROI and scans be in the same space
        if any(any( (roi2fnirtInput - eye(4))>1e-6))
          mrWarnDlg(['(fslApplyWarpROI) Roi ' thisWarpedRoi.name ' is not in the FNIRT input base space, converting ...']);
          thisWarpedRoi.xform = fnirtInputXform;
          thisWarpedRoi.voxelSize = fnirtInputVoxelSize;
          if baseNum>0
            thisRoi.coords = getROICoordinates(thisView,roiList(iRoi),0,[],'baseNum',baseNum);
            thisWarpedRoi.sformCode = viewGet(thisView,'baseSformCode',baseNum);
            thisWarpedRoi.vol2mag = viewGet(thisView,'baseVol2mag',baseNum);
            thisWarpedRoi.vol2tal = viewGet(thisView,'baseVol2tal',baseNum);
          else
            thisRoi.coords = getROICoordinates(thisView,roiList(iRoi));
            thisWarpedRoi.sformCode = viewGet(thisView,'sformCode',scanNum);
            thisWarpedRoi.vol2mag = viewGet(thisView,'scan2mag',scanNum);
            thisWarpedRoi.vol2tal = viewGet(thisView,'scan2tal',scanNum);
          end
        end

        if isempty(thisRoi.coords)
          mrwarnDlg(['Roi ' thisRoi.name ' is empty after conversion, skipping ...']);
          fprintf('\n');
        else
          warpIndices = sub2ind(fnirtInputDims(1:3),thisRoi.coords(1,:),thisRoi.coords(2,:),thisRoi.coords(3,:));
          roiVolume = zeros(fnirtInputDims(1:3));
          roiVolume(warpIndices) = iRoi;
          warpedRoiVolume = fslApplyWarp(roiVolume, [warpCoefPathname warpCoefFilename], tempFilename, fnirtInputHdr, 'nearest');
          warpedCoordsIndices = find(warpedRoiVolume==iRoi);
          [x,y,z] = ind2sub(fnirtInputDims(1:3),warpedCoordsIndices);
          thisWarpedRoi.coords = [x y z ones(length(warpedCoordsIndices),1)]';
        end
      end
    else
      %%%%%%%%%%%%%%%%%%%%%%%% using warp field, but that doesn't take the ROI volume into account
      %it could certainly be modified to do so, but would probably be much slower than using applyWarp, 
      %since taking the volume into account is equivalent to making the ROI continuous, and that's what the other method naturally deals with
      % Also, this has not been re-tested after the xformROIcoords call was moved from above to here
      if any(any( (thisRoi.xform - fnirtInputXform)>1e-6))
        mrWarnDlg(['(fslApplyWarpROI) ROI ' thisRoi.name ' is not in the requested base space, converting ...']);
        scalingFactor = ceil(fnirtInputVoxelSize./thisRoi.voxelSize);
        %transformation matrix between roi space and fnirt input volume space (base)
        roi2fnirtInput = diag([scalingFactor 1])*(fnirtInputXform\thisRoi.xform);
        thisRoi.coords = xformROIcoords(thisRoi.coords,roi2fnirtInput,thisRoi.voxelSize,fnirtInputVoxelSize./scalingFactor);
        fnirtInputDims = scalingFactor.*fnirtInputDims;
        fnirtInputHdr.pixdim(2:4) = fnirtInputHdr.pixdim(2:4)./scalingFactor';
        fnirtInputHdr.sform44 = fnirtInputHdr.sform44*diag([1./scalingFactor 1]);
        fnirtInputHdr.qform44 = fnirtInputHdr.qform44*diag([1./scalingFactor 1]);
      else
        scalingFactor = [1 1 1];
      end

      thisWarpedRoi.coords = fslApplyWarpCoords([thisRoi.coords;ones(1,size(thisRoi.coords,2))], fnirtInputVoxelSize./scalingFactor, [.5 .5 .5], [warpCoefPathname warpCoefFilename], tempFilename, fnirtInputHdr);

      if any(any( (thisRoi.xform - fnirtInputXform)>1e-6))
        %convert back to roi space
        thisWarpedRoi.coords = xformROIcoords(thisWarpedRoi.coords,inv(roi2fnirtInput),fnirtInputVoxelSize./scalingFactor,thisRoi.voxelSize);
      end
    end
    thisView = viewSet(thisView,'newROI',thisWarpedRoi);
    disp(['(fslApplyWarpROI) Warped ROI ' thisRoi.name]);
    fprintf('\n');
    needToRefresh = 1;
  end
end



if needToRefresh
 % select the same current roi
 thisView = viewSet(thisView,'curROI',viewGet(thisView,'roinum',currentROIName));
 refreshMLRDisplay(viewGet(thisView,'viewNum'));
end
   
