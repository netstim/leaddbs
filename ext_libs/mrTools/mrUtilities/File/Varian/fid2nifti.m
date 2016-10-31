% fid2nifti.m
%
%        $Id$ 
%      usage: fid2nifti(fidname,<verbose>,<movepro>,<receiverNum>)
%         by: justin gardner
%       date: 04/30/09
%    purpose: convert a Varian fid into a valid Nifti file with slice orientation information preserved in header.
%             When called with a fidname, it will save a nifti hdr/img pair out with the same name.
%          
%       e.g.: fid2nifti('2danat.fid')
%
%             Can be used with wildcards to process multiple files:
% 
%             fid2nfiti 2d*.fid
% 
%             If you want to just load (not save a converted file)
%
%             [d h] = fid2nifti('fidname.fid');
%
%             You can get verbose information by doing
%
%             [d h] = fid2nifti('fidname.fid',1);
%
%             You can move the pro
%
%             [d h] = fid2nifti('fidname.fid','movepro=4');
% 
%             You can also extract a single receiver in multi-receiver data
%
%             [d h] = fid2nifti('fidname.fid','receiverNum=3');
%
%             or return data with all receivers (default is sum-of-squares)
%     
%             [d h] = fid2nifti('fidname.fid','receiverNum=inf');
%
%             and extract a portion of the image (with the correct
%             header info to match):
%
%             fid2nifti('fidname.fid','xMin=10','xMax=100','yMin=15','yMax=120','sMin=3','sMax=40');
%
%             And specify the outputName to save to
% 
%             fid2nifti('fidname.fid','outputName.hdr');
% 
%             or if you have more arguments, specify output name as:
%
%             fid2nifti('fidname.fid','outputName=outputName.hdr','xMin=10','xMax=100','yMin=15','yMax=120','sMin=3','sMax=40');
%
function [outdata outhdr] = fid2nifti(fidname,varargin)

outdata = [];
outhdr = [];

% check arguments
if nargin < 1
  help fid2nifti
  return
end

% if the first argument is a numeric it means to set verbose setting
verbose = 0;
if (length(varargin) >= 1) && isscalar(varargin{1})
  verbose = varargin{1};
  temp = {};
  for i = 2:length(varargin)
    temp{i-1} = varargin{i};
  end
  varargin = temp;
end

% if there are any equals signs then it is an argument, so split
% the varargin into an argument list and a fidlist
altFidnames = {};argsList = {};
for i = 1:length(varargin)
  if ~isstr(varargin{i}) || isempty(strfind(varargin{i},'.fid'))
    argsList{end+1} = varargin{i};
  else
    altFidnames{end+1} = varargin{i};
  end
end

% parse arguments
movepro=0;receiverNum=[];loadArgs = [];movepss=0;keepref = 0;
xMin=1;xMax=inf;yMin=1;yMax=inf;sMin=1;sMax=inf;outputName=[];kspace=0;rescale=255;
validArgs = {'movepro','movepss','receiverNum','xMin','xMax','yMin','yMax','sMin','sMax','outputName','kspace','rescale','keepref'};
getArgs(argsList,{validArgs{:} 'loadArgs'});
% now evaluate the load args (these are passed from mlrImageLoad
getArgs(loadArgs,validArgs);

% if this has no path, then check for search pattern
if ~iscell(fidname) && strcmp(getLastDir(fidname),fidname)
  % use regexp to serach for matching filenames
  dirlist = dir;
  % make sure fid ends in .fid
  fidname = setext(fidname,'fid',0);

  searchstr = sprintf('^%s$',fixBadChars(fidname,{'*','.*'}));
  fidname = {};
  for i = 1:length(dirlist)
    if regexp(dirlist(i).name,searchstr)
      fidname{end+1} = dirlist(i).name;
    end
  end
end

% cat with other arguments and make sure we have a cell array
if ~isempty(altFidnames)
  fidnames = cellcat(fidname,altFidnames);
else
  fidnames = cellArray(fidname);
end

% make sure we still have something to do
if isempty(fidnames)
  disp(sprintf('(fid2nifti) No matching files'));
  return
end

for i = 1:length(fidnames)
  % get this fidname and make sure it ends in fid
  fidname = fidnames{i};
  fidname = setext(fidname,'fid',0);
  
  if ~isdir(fidname)
    disp(sprintf('(fid2nifti) WARNING: Could not find file %s',fidname));
    continue
  end

  % see if we have to shift pss
  [xform info] = fid2xform(fidname,0);
  % get the fid
  fid = getfid(fidname,'verbose',verbose,'movepro',movepro,'kspace',kspace,'movepss',movepss+info.movepss);
  if isempty(fid.data)
    disp(sprintf('(fid2nifti) WARNING file %s could not be read',fidname));
    continue
  end

  % check if we are being asked to return only a single receiver
  doMergeCoils = true;
  if ~isempty(receiverNum)
    if isinf(receiverNum)
      doMergeCoils = false;
    elseif (receiverNum < 1) || (receiverNum > size(fid.data,5))
      disp(sprintf('(fid2nifti) Data has %s receivers, cannot extract receiver: %i',size(fid.data,5),receiverNum));
      disp(sprintf('            Ignoring receiverNum setting and continuing processing'));
    else
      % extract that receiver
      fid.data = fid.data(:,:,:,:,receiverNum);
      fid.dim = size(fid.data);
    end
  end
  
  % check to see if we have to merge coils
  if (size(fid.data,5) > 1) && ~kspace && doMergeCoils
    numReceivers = size(fid.data,5);
    % see if we have to merge coils
    if numReceivers > 1
      if numReceivers ~= fid.dim(end)
	disp(sprintf('(fid2nifti) Num receivers (%i) does not match data dim (%i)',numReceivers,fid.dim(end)));
      end
      % display what we are doing
      if verbose,disppercent(-inf,sprintf('(fid2nifti) Taking sum of squares of %i coils',numReceivers));end
      % merge the coils
      for volNum = 1:size(fid.data,4)
	sumOfSquares = zeros(fid.dim(1:3));
	for receiverNum = 1:numReceivers
	  sumOfSquares = sumOfSquares+fid.data(:,:,:,volNum,receiverNum).^2;
	end
	data(:,:,:,volNum) = sqrt(sumOfSquares);
	if verbose,disppercent(volNum/size(fid.data,4));end
      end
      if verbose,disppercent(inf);end
      fid.data = data;
      fid.dim = size(data);
    end
  end
  
  % create a header
  hdr = fid2niftihdr(fidname,verbose,'movepro',movepro,'movepss',movepss);
  
  % reorder slices if necessary. note that this is here for fixing interleaved slices,
  % but also reorders slices for 3d images (which go in descending rather than ascending
  % order of pss). This appears not necessary for vnmrs console, so skip
  if ~(strcmp(fid.info.console,'vnmrs') && fid.info.acq3d)
    [pss sliceIndex] = sort(fid.info.procpar.pss);
    if ~isequal(sliceIndex,1:size(fid.data,3))
      fid.data = fid.data(:,:,sliceIndex,:);
    end
  end
  
  % remove any reference volume
  if fid.info.nRefVolumes && ~keepref
    if verbose
      disp(sprintf('(fid2nifti) Removing %i reference volumes',fid.info.nRefVolumes));
    end
    fid.data = fid.data(:,:,:,fid.info.nRefVolumes+1:end,:);
  end
  
  % adjust dimensions if asked for
  [fid.data hdr] = adjustDims(fid.data,hdr,xMin,xMax,yMin,yMax,sMin,sMax);

  % check the dimensions of the data versus the dimensions in the header
  if ~isequal([size(fid.data,1) size(fid.data,2) size(fid.data,3)],hdr.dim(2:4)')
    disp(sprintf('(fid2nifti) Header info from procpar does not match size %s with data read %s',mynum2str(hdr.dim(2:4)'),mynum2str(size(fid.data))));
  end

  % rescale data
  if ~isempty(rescale) && ~isequal(rescale,0)
    if ~kspace
      minData = min(fid.data(:));
      maxData = max(fid.data(:));
      fid.data = rescale*(fid.data-minData)/(maxData-minData);
    end
  end
  
  % write the file, but only if we aren't taking an output argument
  if nargout == 0
    % make a filename 
    if isempty(outputName)
      outputName = setext(fixBadChars(stripext(fidname),{'.','_'}),mrGetPref('niftiFileExtension'));
    else
      outputName = setext(outputName,mrGetPref('niftiFileExtension'));
    end
    disp(sprintf('(fid2nifti) Converting %s to %s',fidname,outputName));
    cbiWriteNifti(outputName,fid.data,hdr);
    outputName = [];
  else
    outhdr{i} = hdr;
    outdata{i} = fid.data;
  end
end

% return single header if that is all that is asked for
if length(outhdr) == 1,
  outdata = outdata{1};
  outhdr = outhdr{1};
end


% cellcat.m
%
%      usage: cellcat()
%         by: justin gardner
%       date: 06/24/07
%    purpose: cat two cell arrays. this preserves order and keeps duplicates (unlike union)
%
function c = cellcat(c1,c2)

% check arguments
if ~any(nargin == [2])
  c={};
  help cellcat
  return
end

% make sure that c1 and c2 are cell arrays
c1 = cellArray(c1);
c2 = cellArray(c2);

c = c1;
for i = 1:length(c2)
  c{end+1} = c2{i};
end

%%%%%%%%%%%%%%%%%%%%
%%   adjustDims   %%
%%%%%%%%%%%%%%%%%%%%
function [data hdr] = adjustDims(data,hdr,xMin,xMax,yMin,yMax,sMin,sMax)

% if nothing to do, just return
if (xMin == 1) && (xMax == inf) && (yMin == 1) && (yMax == inf) && (sMin == 1) && (sMax == inf)
  return
end

% get current dimensions
dims = size(data);

% make the dimensions valid
xMin = round(max(1,xMin));
xMax = round(min(xMax,dims(1)));
yMin = round(max(1,yMin));
yMax = round(min(yMax,dims(2)));
sMin = round(max(1,sMin));
sMax = round(min(sMax,dims(3)));

% adjust data size
data = data(xMin:xMax,yMin:yMax,sMin:sMax,:);

% find out how much the new xMin, yMin, sMin have translated the image
t = hdr.qform44 * [xMin yMin sMin 1]' - hdr.qform44 * [1 1 1 1]';
t = [zeros(3,3) t(1:3,1); 0 0 0 0];

% and compute the new qform
hdr = cbiSetNiftiQform(hdr,t+hdr.qform44);

% reset the dims
dims = size(data);
hdr.dim(2:length(dims)+1) = dims;

