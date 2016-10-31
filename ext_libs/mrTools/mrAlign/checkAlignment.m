% function [inplanes,volume,mosaic] = checkAlignment(correction)
%
%        $Id$
%
% Interpolates the inplanes corresponding to the estimated rotation (rot)
% and translation (trans) and displays slice by slice the original
% inplanes, the interpolated inplanes and a mosaic of both.
%
%     inputs: correction: 0= no correction, 1=intensity/contrast correction, 2=IC correction and reverse contrast
%
% Oscar Nestares - 5/99
% djh 9/2003, modified into a function and cleaned up various things
% djh 1/2007, updated to mrAlign-4.5
% jb 10/2010  - only re-computes volumes if not computed before + implemented volume reduction for fine resolution
%             - put various checks and call to checkAlignmentGUI back into mrAlignGUI
%             - function returns (reduced) inplanes, interpolated volume and mosaic
%             - removed input parameters inplanes, volume and replace by global variable ALIGN
%             - added correction=2: IC correction AND reverse contrast
%
% Julien Besle 10/2010, adapted from checkAlignment by Oscar Nestares
%

function [inplanes,volume,mosaic] = checkAlignment(correction)

global ALIGN

hmsgbox = mrMsgBox('Please wait: interpolating slices, correcting intensity & contrast...');
%reduce volume if necessary
[volume, volumeReductionMatrix] = mrReduceVolume(ALIGN.volume,ALIGN.volumeVoxelSize, 0, 'Volume');
%perform correction on the original inplanes
if correction
   % intensity & contrast correction
   inplanes = intensityContrastCorrection(ALIGN.inplanes, ALIGN.crop);
else
   inplanes = ALIGN.inplanes;
end
%reduce inplanes if necessary
if correction==2
   minInplaneVoxelSize = ALIGN.minEPIFineVoxelSize;
else
   minInplaneVoxelSize = ALIGN.minFineVoxelSize;
end

[inplanes, inplanesReductionMatrix] = mrReduceVolume(inplanes,ALIGN.inplaneVoxelSize, minInplaneVoxelSize, 'Inplanes');
%    ALIGN.xformICCorrection = ALIGN.xform * ALIGN.guiXform;
xform =  volumeReductionMatrix \ ALIGN.guiXform * ALIGN.xform * inplanesReductionMatrix;

% interpolate volume to inplane
volume = interpVolume(volume, xform, size(inplanes));

if correction 
   % volume intensity & contrast correction 
   crop = ALIGN.crop;
   if ~isempty(ALIGN.crop)
      reductionFactor = diag(inplanesReductionMatrix);
      crop = floor(ALIGN.crop./repmat(reductionFactor(1:3)',2,1));
   end
   volume = intensityContrastCorrection(volume, crop);

   if correction==2 %this means we want to reverse contrast
       limit = 4; % threshold used by intensityContrastCorrection
       inplanes = - inplanes * (pi/limit);
       volume = volume * (pi/limit);
   end

   %only keep the computed volumes if correction
   ALIGN.xformICCorrection = ALIGN.guiXform * ALIGN.xform;
   ALIGN.correctedInplanes = inplanes;
   ALIGN.correctedVolume = volume;
   ALIGN.reversedContrast = correction-1;
end

mrCloseDlg(hmsgbox);
   % make mosaic
mosaic = imageMosaic(volume,inplanes);

