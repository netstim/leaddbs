% makescm.m
%
%      usage: makescm(d,<hdrlen>,<applyFiltering>)
%         by: justin gardner
%       date: 07/28/04
%       e.g.: makescm(d)
%    purpose: makes a stimulation convolution matrix
%             for data series. this correctly handles
%             run boundaries. it uses the stimulus volumes
%             found in d.stimvol. if hdrlen is not passed in
%             then it uses d.hdrlen (if it exists). hdrlen is in volumes.
%             applyFiltering defaults to 0, but if it
%             is set will apply the hipass filtering or
%             projection in d.concatInfo to the columns
%             of the scm. d should have the following fields:
%             {dim,stimvol,<concatInfo>,<hdrlen>}
%
%             Can also be used by passing in view and stimvol
%             v = newView;
%             v = viewSet(v,'curGroup','Concatenation'):
%             v = viewSet(v,'curScan',1);
%             stimvol= getStimvol(v,'orientation');
%             scm = makescm(v,20,0,stimvol);
%
function d = makescm(d,hdrlen,applyFiltering,stimvol)

% check arguments
if ~any(nargin == [1 2 3 4])
  help makescm
  return
end

% deal with calling by view instead of d structure
returnOnlySCM = 0;
if isview(d)
  v = d;d = [];
  d.concatInfo = viewGet(v,'concatInfo');
  d.dim = viewGet(v,'scanDims');
  % get the length of scan
  d.dim(end+1) = viewGet(v,'nFrames');
  d.concatInfo = viewGet(v,'concatInfo');
  d.nframes = viewGet(v,'nFrames');

  if isempty(stimvol)
    disp(sprintf('(makescm) Must pass in stimvol'));
    return
  else
    d.stimvol = stimvol;
  end
  returnOnlySCM = 1;
end
  

% set hdrlen
if ieNotDefined('hdrlen')
  if (isfield(d,'hdrlen'))
    hdrlen = d.hdrlen;
  else
    hdrlen = 25;
  end
end

if isfield(d,'nFrames') && (hdrlen > d.nFrames)
    mrErrorDlg(sprintf('(makescm) Your hdrlen in volumes (%i) is greater than the length of your scan (%i), perhaps framePeriod is set incorrectly?',hdrlen,d.nFrames));
end

% default to not apply filtering that concatTSeries did
if ieNotDefined('applyFiltering'), applyFiltering = 0;end

% if we have only a single run then we set
% the runTransitions for that single run
if ~isfield(d,'concatInfo') || isempty(d.concatInfo)
  runTransition = [1 d.dim(end)];
else
  runTransition = d.concatInfo.runTransition;
end

% go through each run of the experiment
%allscm = [];
allscm = zeros(runTransition(end,2),length(d.stimvol)*hdrlen);
for runnum = 1:size(runTransition,1)
  % default values
%  scm = [];
  scm = zeros(runTransition(runnum,2)-runTransition(runnum,1)+1,length(d.stimvol)*hdrlen);;
  % make stimulus convolution matrix
  for stimnum = 1:length(d.stimvol)
    % make an array containing the stimulus times
    stimarray = zeros(1,runTransition(runnum,2)-runTransition(runnum,1)+1);
    % only use stimvols that are within this runs volume numbers
    stimarray(d.stimvol{stimnum}(find((d.stimvol{stimnum}>=runTransition(runnum,1)) & (d.stimvol{stimnum}<=runTransition(runnum,2))))-runTransition(runnum,1)+1) = 1;
    % stack stimcmatrices horizontally
    m = stimconv(stimarray,hdrlen);
    % apply the same filter as original data
    if applyFiltering
      % check for what filtering was done
      if isfield(d,'concatInfo') 
	% apply hipass filter
	if isfield(d.concatInfo,'hipassfilter') && ~isempty(d.concatInfo.hipassfilter{runnum})
	  m = real(ifft(fft(m) .* repmat(d.concatInfo.hipassfilter{runnum}', 1, size(m,2)) ));
	end
	% project out the mean vector
	if isfield(d.concatInfo,'projection') && ~isempty(d.concatInfo.projection{runnum})
	  projectionWeight = d.concatInfo.projection{runnum}.sourceMeanVector * m;
	  m = m - d.concatInfo.projection{runnum}.sourceMeanVector'*projectionWeight;
	end
	% now remove mean
	m = m-repmat(mean(m,1),size(m,1),1);
      end
    end
    scm(:,((stimnum-1)*hdrlen+1):stimnum*hdrlen) = m;
%    scm = [scm, m];
  end
  % stack this run's stimcmatrix on to the last one
  allscm(runTransition(runnum,1):runTransition(runnum,2),:) = scm;
%  allscm = [allscm; scm];
end

% set values
d.nhdr = length(d.stimvol);
d.scm = allscm;
d.hdrlen = hdrlen;
d.volumes = 1:d.dim(end);

if returnOnlySCM
  d = d.scm;
end

