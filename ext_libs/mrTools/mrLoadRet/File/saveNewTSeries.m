function [view,filename] = saveNewTSeries(view,tseries,scanParams,hdr,makeLink)
%
% view = saveNewTSeries(view,tseries,[scanParams],[hdr],[makeLink])
%
% tseries: tseries array (x,y,z,t) or a filename of a nifti file
%
% scanParams: structure with the following fields: fileName, description,
% junkFrames, nFrames. The other scanParams fields are set automatically
% from the nifti header. 
% Default:
%    scanParams.fileName = 'tseries-mmddyy-hhmmss.img'; 
%       where mmddyy = date and hhmmss = time
%    scanParams.junkFrames = 0;
%    scanParams.nFrames = size(tseries,4);
%    scanParams.description = '';
%
% Chooses between .img and .nii based on 'niftiFileExtension' preference.
%
% hdr: template for nifti header. The header is always passed through
% cbiCreateNiftiHeader to ensure consistency with the data. Default: [];
%
% makeLink: Optional, if tseries is passed as a filename then setting
% this to 1 will cause the tseries to be linked rather than copied
% defaults to 0
%
% $Id$
%

% check arguments
if ~any(nargin == [1:5])
  help saveNewTSeries
  return
end



if ieNotDefined('scanParams')
  scanParams.fileName = [];
  scanParams.junkFrames = 0;
  scanParams.description = '';
end
if ieNotDefined('hdr')
  hdr = [];
end
if ieNotDefined('makeLink'),makeLink=0;end

if isempty(scanParams.fileName)
  if ischar(tseries)
    extension = ['.',getext(tseries)];
  else
    extension = mrGetPref('niftiFileExtension');
  end
  scanParams.fileName = ['tseries-',datestr(now,'yymmdd'),'-',datestr(now,'HHMMSS'),extension];
end
filename = scanParams.fileName;

% Save tseries 
tseriesdir = viewGet(view,'tseriesdir');
path = fullfile(tseriesdir,scanParams.fileName);

% see if the tseries is actually a string in which case we should
% copy the nifti file.
if ischar(tseries)
  success = copyNiftiFile(tseries,path,makeLink);
  if ~success,return,end

  % get the number of frames
  hdr = mlrImageReadNiftiHeader(tseries);
  nFrames = hdr.dim(5);
else
  [byteswritten,hdr] = cbiWriteNifti(path,single(tseries),hdr,'float32');
  nFrames = size(tseries,4);
end

% set nFrames
if ~isfield(scanParams,'nFrames')
  scanParams.nFrames = nFrames;
end  
% Add new scan (note that the tseries file must already exist)
view = viewSet(view,'newScan',scanParams);
