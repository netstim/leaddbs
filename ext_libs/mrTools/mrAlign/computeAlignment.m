function outputXform = computeAlignment(resolution, handles)
% xform = computeAlignment()
%
% Oscar Nestares - 5/99
%
%        $Id$
%
% ON 5/00     - NEW PARAMETER NzLIMIT that controls if slices are replicated or
%               not in the edges. I've put a high value so that almost ALWAYS the
%               slices are replicated before doing the alignment.
%             - Added an alarm when maximum number of iterations is reached, with
%               the possibility of continuing the iterations.
% SL 7/25/02 -  Created regVolInp4 for use in mrAlign4 (based entirely on regVolInp)
% SL 7/31/02 -  Changed instances of A*inv(B), where A and B are both matrices,
%               to be A/B.
% SL 8/01/02 -  Added coarseFlag and fineFlag options
% djh 9/2003 -  Cleaned up various things, fixed a major bug (although I'm
%               entirely sure what it was).
% djh 5/2004 -  Added crop.
% djh 1/2007 -  Update to mrAlign-4.5
% djh 7/2007 -  Added contrastReverse param
% jb 10/2010 -  Added ignoreZeroVoxels parameter: changes zeros to NaNs after intensity correction so that these voxels do not influence motion estimation (at least not too much)
%            -  Included volume reduction, parametrized by minVolumeVoxelSize and minInplanesVoxelSize: fixed minimum voxel size in millimeters
%            -  removed all parameters, replaced by the use of global structure ALIGN
%            -  added parameter handles: handles to the figure to update the GUI on the fly
%            -  added parameter resolution: 'fine' or 'coarse'

outputXform = [];

mindisp = 0.1;  % displacement used to end the iterations
maxDispLimit = 50; %displacement used to detect abnormal behaviour of the estimation
CB = 2.5;       % cutoff parameter for robust estimation
SC = 2;         % scaling parameter for robust estimation
nzLIMIT = 16;   % if the number of inplanes is less that this, then
                % it puts a border by replicating the first and last
                % inplanes two times.
permuteM = [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];  % For swapping x and y

global ALIGN;

if ALIGN.currentlyComputingAlignment
  mrWarnDlg('(computeAlignment) An alignment is already being computed')
  return;
end

%Perform various checks
if isempty(ALIGN.volume) || isempty(ALIGN.inplanes)
	mrWarnDlg('(computeAlignment) Load source and destination volumes before computing alignment');
	return
end
if isempty(ALIGN.xform) 
	mrWarnDlg('(computeAlignment) Initialize aligment or load a previously saved alignment before computing.');
	return
end

if ieNotDefined('resolution'),resolution = 'fine' ;end
switch(resolution)

   case 'fine'
      minVolumeVoxelSize = 0;
      switch(ALIGN.reverseContrastAlignment)
         case 0
            minInplanesVoxelSize = ALIGN.minFineVoxelSize;
         case 1
            minInplanesVoxelSize = ALIGN.minEPIFineVoxelSize;
      end
      
   case 'coarse'
      switch(ALIGN.reverseContrastAlignment)
         case 0
            minVolumeVoxelSize = ALIGN.minCoarseVoxelSize;
            minInplanesVoxelSize = ALIGN.minCoarseVoxelSize;
         case 1
            minVolumeVoxelSize = ALIGN.minEPICoarseVoxelSize;
            minInplanesVoxelSize = ALIGN.minEPICoarseVoxelSize;
      end
end


crop = ALIGN.crop;
%ALIGN.inplanes = ALIGN.inplanes;
%get the current xform
oldXform = ALIGN.xform;
[oldGuiTranslation oldGuiRotation] = getAlignGUI(handles);
xform = ALIGN.guiXform * ALIGN.xform;
if ALIGN.rigidBodyAlignment && ~ALIGN.xformIsRigidBody   %if rigid-body alignment asked but current xform is non-rigid, ask if continue
   answer = questdlg(...
      'Non-rigid body has previously been performed and current xform might include non-rigid body transformations. Are you sure you want to continue ?',...
      'Rigid body after non-rigid body','Yes','No','No');
   if strcmp(answer,'No')
      return
   end
end  

%reduce source volume, if needed
[volume, volumeReductionMatrix] = mrReduceVolume(ALIGN.volume,ALIGN.volumeVoxelSize, minVolumeVoxelSize, 'Volume');

%if ignore zero voxels
if ALIGN.ignoreZeroVoxels     %mask has to be computed before intensity/contrast correction
  inplanesMask = ALIGN.inplanes;
  inplanesMask(isnan(inplanesMask)) = 0;
  inplanesMask = logical(inplanesMask);
else
  mask = [];
end

% correct intensity & contrast of inplane. Results are better if done on non-reduced inplanes
wbh = mrMsgBox('Please wait: interpolating slices, correcting intensity & contrast...');
inpIC = intensityContrastCorrection(ALIGN.inplanes, crop);
mrCloseDlg(wbh);

%reduce destination volume, if needed     
[inpIC, inplanesReductionMatrix] = mrReduceVolume(inpIC,ALIGN.inplaneVoxelSize, minInplanesVoxelSize, 'Inplanes');
xform =  volumeReductionMatrix \ xform * inplanesReductionMatrix;
if ~isempty(crop)
   reductionFactor = diag(inplanesReductionMatrix);
   crop = floor(crop./repmat(reductionFactor(1:3)',2,1));
end

if ALIGN.ignoreZeroVoxels %we need to reduce the mask
   inplanesMask = mrReduceVolume(inplanesMask+0,ALIGN.inplaneVoxelSize, minInplanesVoxelSize, 'Mask');
   inplanesMask = inplanesMask>.8;  %this value could be increased to be more conservative 
                                    %and exclude voxels that have been modified by IC correction
end
                                                
if ALIGN.reverseContrastAlignment
    limit = 4; % threshold used by intensityContrastCorrection
    inpIC = inpIC * (pi/limit);
    inpIC = -inpIC;
end

% If the number of slices is too small, repeat the first and last slice
% to avoid running out of data (the derivative computation discards the
% borders in z, 2 slices at the begining and 2 more at the end)
inplanesSize = size(inpIC);
if inplanesSize(3) < nzLIMIT
	handleMsg = mrMsgBox(['Only ',num2str(size(inpIC,3)),' slices. Padding by replicating first and last slice.']);
	inpIC = cat(3,inpIC(:,:,1),inpIC(:,:,1),inpIC,...
		inpIC(:,:,end),inpIC(:,:,end));
   if ALIGN.ignoreZeroVoxels
      inplanesMask = cat(3,inplanesMask(:,:,1),inplanesMask(:,:,1),inplanesMask,...
         inplanesMask(:,:,end),inplanesMask(:,:,end));
   end
   if ~isempty(crop)
     crop(:,3) = crop(:,3)+2;
   end
   mrCloseDlg(handleMsg);
end

% Loop until the approximate maximum displacement is less than MINDISP, or
% the maximum number of iterations is reached. The maximum displacement is
% calculated aproximately from the sum of terms: the displacement
% corresponding to the rotation of the farthest point in the inplanes, plus
% the norm of the translation.
keepIterating = 'Yes';
maxdisp = mindisp + 1;
maxIter = ALIGN.NIter;
ALIGN.currentlyComputingAlignment=1;
ALIGN.stopComputingAlignment = 0;
saveNewXform = 'Yes';

while strcmp(keepIterating, 'Yes') && ~ALIGN.stopComputingAlignment
   niter = 0;
   wbh = mrWaitBar(0,'Computing alignment...');
   while ( ((maxdisp > mindisp) && (niter < maxIter)) || (niter==0)) && ~ALIGN.stopComputingAlignment 

     % updating number of iterations and message in progress bar
     niter = niter+1;
     mrWaitBar(niter/maxIter,wbh);

     % resample volume according to initial xform
     volInterp = interpVolume(volume, xform, inplanesSize);

     % Pad additional slices, if needed
     if inplanesSize(3) < nzLIMIT
         volInterp = cat(3,volInterp(:,:,1),volInterp(:,:,1),volInterp,...
             volInterp(:,:,end),volInterp(:,:,end));
     end
     if ALIGN.ignoreZeroVoxels %if ignore zero voxels
        mask = volInterp;
        mask(isnan(mask)) = 0;
        mask = logical(mask) & inplanesMask;
     end
     % correct intensity & contrast of interpolated volume
     try
         volIC = intensityContrastCorrection(volInterp, crop);
     catch exception
        mrWarnDlg(['A problem has occurred during intensity/contrast correction. ' exception.message '. Aborting ...']);
        saveNewXform = 'No';
        break;
     end

     if ALIGN.reverseContrastAlignment
         volIC = volIC * (pi/limit);
     end
     % *** For debugging ***
     % xform
     % figure(1); imagesc(inpIC(:,:,9)); colormap(gray); axis image;
     % figure(2); imagesc(volIC(:,:,9)); colormap(gray); axis image;

     % motion estimation (no multiresolution, no iterative)
     M = estMotion3(inpIC, volIC,...
         ALIGN.rigidBodyAlignment,...	      % rotFlag
         ALIGN.robustAlignment,...% robustFlag 
         ALIGN.reverseContrastAlignment,... % treat image intensities as phase valued
         crop,...      % crop
         CB,...        % cutoff parameter for robust estimation
         SC,...        % scale parameter for robust estimation
         mask);        % mask for voxels that should not be taken into account in the minimization 
     % estMotion3 uses [y,x,z] convention so we need to swap x and y
     % before passing coordinates through M and then swap them back
     % before passing them through xform. Yuch. But not worth rewriting
     % estMotion3.
     M = permuteM * M * permuteM;
     xform = xform * M;
     maxdisp = norm(inplanesSize' - M(1:3,1:3)*inplanesSize') + norm(M(1:3,4));

     if ALIGN.displayAlignmentSteps%refresh display
        ALIGN.xform = ALIGN.guiXform \volumeReductionMatrix* xform / inplanesReductionMatrix;
        refreshAlignDisplay(handles);
     end

     if maxdisp>maxDispLimit
        question = ['The estimated displacement seems abnormally large (' num2str(maxdisp) 'mm). '];
        question = [question 'If problem persists, try Ignoring zero voxels or switching off Robust alignment. '];
        question = [question 'Do you want to continue iterating? '];
        keepIterating = questdlg(question, 'Continue Iterating?', 'Yes','No','No');
        if strcmp(keepIterating,'No')
            saveNewXform = 'No';
            break;
        end
     end
   end
   if (niter == maxIter)
     mrCloseDlg(wbh);
     % ask if we should continue iterating
     quest = {['Maximum number of iterations (',num2str(maxIter),') reached.'], 'Continue Iterating?'};
     keepIterating = questdlg(quest, 'WARNING', 'Yes', 'No', 'No');
   else
     if strcmp(saveNewXform,'No') && isstruct(wbh)
        wbh = rmfield(wbh,'disppercent');
     end
     mrCloseDlg(wbh);
     keepIterating = 'No';
   end
    
end
ALIGN.currentlyComputingAlignment =0;


if ALIGN.stopComputingAlignment
   ALIGN.stopComputingAlignment = 0;
   saveNewXform = questdlg('Do you want to keep the newly computed xform ?', 'Keep Computed Xform ?','Yes','No','No');
end
   
if strcmp(saveNewXform,'Yes')
   ALIGN.oldXform = oldXform;
   ALIGN.oldGuiTranslation = oldGuiTranslation;
   ALIGN.oldGuiRotation = oldGuiRotation;
   ALIGN.xform = volumeReductionMatrix*xform/inplanesReductionMatrix;

   % resample and IC correct volume according to final xform, for later use (in checkAlignment)
   wbh = mrMsgBox('Please wait: interpolating slices, correcting intensity & contrast...');
   volInterp = interpVolume(volume, xform, inplanesSize);
   if inplanesSize(3) < nzLIMIT
      volInterp = cat(3,volInterp(:,:,1),volInterp(:,:,1),volInterp,...
          volInterp(:,:,end),volInterp(:,:,end));
   end
   ALIGN.correctedVolume = intensityContrastCorrection(volInterp, crop);
   if ALIGN.reverseContrastAlignment
      ALIGN.correctedVolume = ALIGN.correctedVolume * (pi/limit);
      ALIGN.reversedContrast = 1;
   else
      ALIGN.reversedContrast = 0;
   end
   ALIGN.correctedInplanes = inpIC;
   mrCloseDlg(wbh);
   
   setAlignGUI(handles,'rot',[0 0 0]);
   setAlignGUI(handles,'trans',[0 0 0]);
   ALIGN.guiXform = getGuiXform(handles);
   if ALIGN.rigidBodyAlignment
      ALIGN.xformIsRigidBody =1;
   else
      ALIGN.xformIsRigidBody =0;
   end
   ALIGN.xformICCorrection = ALIGN.xform;
else
   ALIGN.xform = oldXform;
end
outputXform = ALIGN.xform;

% Refresh display
refreshAlignDisplay(handles);




