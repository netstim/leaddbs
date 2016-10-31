% mrUpSample.m
%
%        $Id$	
%      usage: [output, hdr] = mrUpSample(input, varargin)
%         by: eli merriam
%       date: 08/24/07
%    purpose: 
%
% Upsamples the image im by the integer nLevels.  The upsampled image
% is blurred and edges are dealt with by zero-padding.
%
% The blurring is done with the filter kernel specified by filt
% (default = [.125 .5 .75 .5 .125]'), which should be a vector
% (applied separably as a 1D convolution kernel in X, Y, and Z).
%
% Examples:
%
%     input is a nifti file, output is an upsampled nifti file
%            mrUpSample(inputFileName, 'outName=upSamp');
%
%     input is a nifti file, prompts user for output nifti file name
%            mrUpSample(inputFileName);
%
%     input is a matlab array, output is an upsampled array, and possible a nifti header
%           [upSampleData newHeader] = mrUpSample(data, 'hdr', niftiHeader);
%
% Additinoal arguments:
%     nLevels (default = 1)
%     filt (default = sqrt(2)*namedFilter('binom5'))
%
% 99.08.16 RFD wrote it, based on upBlur (and helpful
%				comments from DJH)
% 07.08.19, EPM added upsampling across slices, nifti file handling
 
function [varargout] = mrUpSample(data, varargin)

% check arguments
if ~any(nargin == [1 2 3])
  help mrUpSample
  return
end

if isstr(data)
  % Input is a file
  isfile = 1;
  [data, hdr] = mlrImageReadNifti(data);
else
  % Input is a Matlab data array
  isarray = 1;
end

% evaluate the arguments
eval(evalargs(varargin,0));

% double the resolution
if ieNotDefined('nLevels'), nLevels = 1; end
if nLevels ~= fix(nLevels)
    error('nLevels must be an integer!!!');
end

% default filter for post-upsample convolution
if ieNotDefined('filt'), filt = sqrt(2)*namedFilter('binom5'); end


% use a little recursion to deal with upsample steps > 2
% ATTN need to fix this.
if nLevels > 1
    data = mrUpSample(data, nLevels-1);
end
if (nLevels >= 1)
    if (any(size(data)==1))
        if (size(data,1)==1)
            filt = filt';
        end
        output = upConv(data, filt, 'zero', (size(data)~=1)+1);
    else
        
        [nrows ncols nslices nframes] = size(data);

        disppercent(-inf,'(mrUpSample) upsampling the rows');
        data = reshape(data,[nrows ncols*nslices*nframes]);
        data = upConv(data, filt, 'zero', [2 1]);
        nrows = nrows*2;
        data = reshape(data,[nrows ncols nslices nframes]);
        disppercent(inf);

        disppercent(-inf,'(mrUpSample) upsampling the columns');
        data = permute(data,[2 1 3 4]);
        data = reshape(data,[ncols nrows*nslices*nframes]);
        data = upConv(data, filt, 'zero', [2 1]);
        ncols = ncols*2;
        data = reshape(data,[ncols nrows nslices nframes]);
        data = permute(data,[2 1 3 4]);
        disppercent(inf);
        
        if size(data,3) > 1
          disppercent(-inf,'(mrUpSample) upsampling the slices');
          data = permute(data,[3 2 1 4]);
          data = reshape(data,[nslices ncols*nrows*nframes]);
          data = upConv(data, filt, 'zero', [2 1]);
          nslices = nslices*2;
          data = reshape(data,[nslices ncols nrows nframes]);
          data = permute(data,[3 2 1 4]);
          disppercent(inf);
        end
    end
else
    data = input;
end

% check if input was nifti, if so, adjust
if ~ieNotDefined('hdr')
    scaleFactor = [nLevels nLevels nLevels]*2;
    % set the qform/sform
    hdr = cbiSetNiftiQform(hdr,hdr.qform44*diag([1./scaleFactor 1]));
    hdr = cbiSetNiftiSform(hdr,hdr.sform44*diag([1./scaleFactor 1]));
end

% save, if the input was a file
if ~ieNotDefined('isfile') & ~ieNotDefined('hdr')
    % set the file extension
    niftiFileExtension = mrGetPref('niftiFileExtension');
    if isempty(niftiFileExtension)
        niftiFileExtension = '.img';
    end
    
    % prompt user for output name
    if ieNotDefined('outName')
        outName = putPathStrDialog(pwd,'Specify name of upsampled Nifti file',['*' niftiFileExtension]);
    end
    
    % write the file
    fprintf('Saving %s...\n', outName);
    [byteswritten,hdr] = cbiWriteNifti(outName, data, hdr);
        
else
    % if input was not a file, return the upsampled data
    varargout{1} = data;
    if ~ieNotDefined('hdr')    
        varargout{2} = hdr;
    end
end

fprintf('done\n');
return;


