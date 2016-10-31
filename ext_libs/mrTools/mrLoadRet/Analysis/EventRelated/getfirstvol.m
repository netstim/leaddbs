% getfirstvol.m
%
%      usage: d = getfirstvol(d)
%         by: justin gardner
%       date: 02/08/05
%
function d = getfirstvol(d)

if (nargin ~= 1)
  help getfirstvol;
  return
end

if (~isfield(d,'stimvol'))
  if isfield(d,'channels')
    % calculate the times of the stimuli
    % get the stim times
    stimraw = d.channels(d.stimchannel,:);
    stimraw(stimraw < 0) = 0;
    stimtimes = find([0 (diff(stimraw)~=0)]);
  
    % get the image number
    acqnum = cumsum(d.acq>1);

    % set the beginning acqnum to 1, so that
    % any event that happens before the first
    % acquistion pulse is assumed to happen
    % during the first acquisition pulse.
    acqnum(1:first(find(acqnum == 1))) = 1;

    % sort into stimuli
    nhdr = max(stimraw);
    for i = 1:nhdr
      stimtimescell{i} = stimtimes(stimraw(stimtimes) == i);
      pulselens(i) = i;
      stimvol{i} = acqnum(stimtimescell{i});
    end
    stimtimes = stimtimescell;
  else
    disp(sprintf('(getfirstvol) No stimvol in d structure'));
    return
  end
else
  stimvol = d.stimvol;
end

% now look in the stimvol structure to see where the
% first stimulus volume is.
if (iscell(stimvol))
  d.firstvol = inf;
  if exist('pulselens','var')
    if ((pulselens(1) ~= 1) && (pulselens(1) < 50)),start = 2;,else,start = 1;,end
  else
    start = 1;
  end
  for i = start:length(stimvol)
    if ~isempty(stimvol{i}) && (stimvol{i}(1) < d.firstvol)
      d.firstvol = stimvol{i}(1);
    end
  end
else
  d.firstvol = stimvol(1);
end

