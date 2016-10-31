function[] = mrExport2SR(viewNum, pathstr)
% mrExport2SR.m
%
%      usage: [] = mrExprt2SR(viewNum, pathstr)
%         by: eli merriam
%       date: 03/20/07
%    purpose: exports a MLR overlay to a Nifti file compatible with SurfRelax
%        $Id$	
%
% modified by julien besle 22/01/2010 to speed up things and take getBaseSpaceOverlay.m out

%mrGlobals

% Get view
view = viewGet(viewNum,'view');

% Get values from the GUI
scanNum = viewGet(view,'curscan');
baseNum = viewGet(view,'currentBase');
overlayNum = viewGet(view,'currentOverlay');
overlayData = viewGet(view,'overlayData',scanNum,overlayNum);


%basedims = viewGet(view, 'basedims');

%transform values in base space
[new_overlay_data, new_base_voxel_size] = getBaseSpaceOverlay(view, overlayData, scanNum, baseNum);

if isempty(new_overlay_data)
  return
end
%write nifti file
baseVolume = viewGet(viewNum,'baseVolume');
hdr = baseVolume.hdr;
%hdr.dim = [size(size(new_overlay_data),2); size(new_overlay_data,1); size(new_overlay_data,2); size(new_overlay_data,3); 1; 1; 1; 1];
hdr.bitpix = 32;   
hdr.datatype = 16;
hdr.is_analyze = 1;
hdr.scl_slope = 1;
hdr.endian = 'l';
if any(new_base_voxel_size ~= viewGet(view,'basevoxelsize',baseNum))
   hdr.pixdim = [0 new_base_voxel_size 0 0 0 0]';        % all pix dims must be specified here
   hdr.qform44 = diag([new_base_voxel_size 0]);
   hdr.sform44 = hdr.qform44;
end
   
% set the file extension
niftiFileExtension = mrGetPref('niftiFileExtension');
if isempty(niftiFileExtension)
  niftiFileExtension = '.img';
end

cbiWriteNifti(sprintf('%s%s',stripext(pathstr),niftiFileExtension), new_overlay_data,hdr);

return







