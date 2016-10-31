% getfidk.m
%
%        $Id:$ 
%      usage: getfidk(filename)
%         by: justin gardner
%       date: 08/10/11
%    purpose: returns raw fid data using getfidkraw.c. Reorders lines of k-space appropriately
%             so that images can be fourier transformed
%
function d = getfidk(filename,varargin)

d.data = [];
% check arguments
if nargin < 1
  help getfidk
  return
end

% parse arguments
verbose = [];
getArgs(varargin,{'verbose=0'});

% check filename
if ~isdir(filename)
  filename = setext(filename,'fid');
  if ~isdir(filename)
    disp(sprintf('(getfidk) Could not find fid directory %s',filename));
    return
  end
end

% load the infromation from the procpar using fid2xform
[xform info] = fid2xform(filename);
if isempty(info),return,end

% look for fid file
fidFilename = fullfile(filename,'fid');
if ~isfile(fidFilename)
  disp(sprintf('(getfidk) Could not find fid file: %s', fidFilename));
  return
end

% load the fid file using getfidkraw
d = getfidkraw(fidFilename,verbose);

% take transpose to make these easier to deal with
d.real = d.real';
d.imag = d.imag';

% compute ow any volumes, slices, phase encode lines and receivers we have
numVolumes = info.dim(4)+info.nRefVolumes;
numSlices = info.dim(3);
numPhaseEncodeLines = info.dim(2)/info.accFactor;
numReceivers = info.numReceivers;

% now compute how many lines we need to read the data from
if info.compressedFid
  % for compressedFid each line of data has all phase encode lines
  numLines = numSlices*numVolumes*numReceivers;
  % but, if intlv is 'y', then each block of data contains all slices
  if isequal(info.intlv,'y')
    numLines = numLines/numSlices;
    if verbose,disp(sprintf('(getfidk) intlv is %s: Assuming shot interleaving',info.intlv));end
  end
  % we will also need to know the line order as returned form petable
  if ~isfield(info.procpar,'pelist') || length(info.procpar.pelist) > 1
    lineorder = info.procpar.pelist-min(info.procpar.pelist)+1;
  else
    % pelist variable not set, read petable
    petable = readpetable(info.procpar.petable{1},'verbose',verbose);
    % petable not found, give up
    if isempty(petable)
      disp(sprintf('(getfidk) !!!! Unknown petable !!!!'));
      keyboard
    end
    % get the line ordering from petable
    lineorder = petable.t1';
    lineorder = lineorder-min(lineorder)+1;
  end
else
  numLines = numPhaseEncodeLines*numSlices*numVolumes*numReceivers;
end

% make sure the dimensions match
if d.nblocks ~= numLines
  disp(sprintf('(getfidk) !!! Number of lines of k-space found in fid (%i) does not match expected number %i (%ix%ix%ix%i %i receivers accFactor %i)',d.nblocks,numLines,info.dim(1),info.dim(2),info.dim(3),info.dim(4),info.numReceivers,info.accFactor));
  disp(sprintf('(getfidk) !!! This is usually because epi processing (epibsi) has not been done !!!'));
  d.data = [];
  return
end

if verbose
  disp(sprintf('(getfidk) Numlines: %i (%ix%ix%ix%i %i receivers accFactor %i)',numLines,info.dim(1),info.dim(2),info.dim(3),info.dim(4),info.numReceivers,info.accFactor));
end

% read the data from fid block structure
kNum = 1;clear i;
d.data = nan(info.dim(1),numPhaseEncodeLines,numSlices,numReceivers,numVolumes);
if verbose,disppercent(-inf,'(getfidk) Reordering data');end
if info.compressedFid
  % if intlv is set to y then it means that shots are interleaved - i.e. a shot is taken on
  % each slice and then you come back. Thus each block contains the data from all slices
  if isequal(info.intlv,'y')
    % etl can be computed
    numshots = info.procpar.numshots;
    etl = numPhaseEncodeLines/numshots;
    % but, make sure we have etl field which tells how many lines per shot
    if ~isfield(info.procpar,'etl') || isempty(info.procpar.etl)
      disp(sprintf('(getfidk) Missing etl field which should tell us how many lines of k-space per shot'));
      disp(sprintf('          Using %i based on having %i shots',etl,info.procpar.numshots));
    else
      % make sure they match what we calculated
      if ~isequal(etl,info.procpar.etl)
	disp(sprintf('(getfidk) Computed etl: %i does not match etl in procpar: %i',etl,info.procpar.etl));
      end
      % get etl from procpar setting
      etl = info.procpar.etl;
      % and get numshots from that
      numshots = numPhaseEncodeLines/etl;
      if numshots ~= round(numshots)
	disp(sprintf('(getfidk) Etl of %i does not divide evenly into num of phase encode lines %i',etl,numPhaseEncodeLines));
	numshots = round(numshots);
      end
    end

    subblockSize = info.dim(1)*etl*numSlices;
    for volNum = 1:numVolumes
      for receiverNum = 1:numReceivers
	blockData = double(d.real(kNum,:)) + i*double(d.imag(kNum,:));
	kNum = kNum + 1;
	for shotNum = 1:numshots
	  % compressed fids with intlv have all lines of k-space in one
	  % single block of data note conversion here to double
	  d.data(:,lineorder((shotNum-1)*etl+1:shotNum*etl),:,receiverNum,volNum) = reshape(blockData((shotNum-1)*subblockSize+1:shotNum*subblockSize),info.dim(1),etl,numSlices);
	end
      end
      if verbose,disppercent(calcPercentDone(volNum,numVolumes,receiverNum,numReceivers));end
    end
  else
    % do processing for non itls compressed data
    for volNum = 1:numVolumes
      for sliceNum = 1:numSlices
	for receiverNum = 1:numReceivers
	  % compressed fids have all lines of k-space in one single block of data
	  % note conversion here to double
	  d.data(:,lineorder,sliceNum,receiverNum,volNum) = reshape(double(d.real(kNum,:)) + i*double(d.imag(kNum,:)),numPhaseEncodeLines,info.dim(2));
	  kNum = kNum+1;
	end
      end
      if verbose,disppercent(calcPercentDone(sliceNum,numSlices,receiverNum,numReceivers,volNum,numVolumes));end
    end
  end
else
  % do processing for non-compressed data
  for kLine = 1:numPhaseEncodeLines
    for volNum = 1:numVolumes
      for receiverNum = 1:numReceivers
	for sliceNum = 1:numSlices
	% uncompress fid contains one line of k-space per block
	  % note conversion here to double
	  d.data(:,kLine,sliceNum,receiverNum,volNum) = double(d.real(kNum,:)) + i*double(d.imag(kNum,:));
	  kNum = kNum+1;
	end
      end
    end
    if verbose,disppercent(calcPercentDone(kLine,numPhaseEncodeLines,volNum,numVolumes,receiverNum,numReceivers,sliceNum,numSlices));end
  end
end

% remove original data read from getfidkraw
d = rmfield(d,'real');
d = rmfield(d,'imag');

if verbose,disppercent(inf);end

