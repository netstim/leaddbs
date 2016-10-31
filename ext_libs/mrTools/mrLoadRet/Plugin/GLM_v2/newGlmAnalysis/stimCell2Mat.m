% stimCell2Mat.m
%
%      usage: stimCell2Mat(stimOnsets, stimDurations{iStim}, framePeriod)
%         by: julien besle 
%       date: 05/12/2010
%    purpose: Transforms stimulus onsets/durations in cell array of vectors form into a stimulus matrix
%              $Id$

function [stimMatrix,runTransitions] = stimCell2Mat(stimOnsets, stimDurations,runTransitions)

if ~ieNotDefined('stimDurations')
  %check that stimDurations and stimOnsets are compatible
  if length(stimDurations)~=length(stimOnsets)
    mrErrorDlg(sprintf('(stimCell2Mat) stimulus onsets and durations are not compatible (%d vs %d event types)',length(stimOnsets),length(stimDurations)));
  end
  %apply stimDuration
  for iStim = 1:length(stimOnsets)
    if length(stimDurations{iStim})~=length(stimOnsets{iStim})
      mrErrorDlg(sprintf('(stimCell2Mat) stimulus onsets and durations are not compatible (%d vs %d event types)',length(stimOnsets{iStim}),length(stimDurations{iStim})));
    end
    if ~isempty(stimOnsets{iStim})
      stimOnsets{iStim} = reshape(stimOnsets{iStim},1,numel(stimOnsets{iStim})); %make sure stimOnsets{iStim} and stimDurations{iStim} 
      stimDurations{iStim} = reshape(stimDurations{iStim},1,numel(stimDurations{iStim})); %are row vectors
      maxDuration = max(stimDurations{iStim});
      stimNum = length(stimDurations{iStim});
      stimPresent = (cumsum(ones(maxDuration,stimNum),1) ./ repmat(stimDurations{iStim},maxDuration,1))<=1;
      stimOnsets{iStim} = repmat(stimOnsets{iStim},maxDuration,1);
      stimOnsets{iStim} = stimOnsets{iStim} + cumsum(stimPresent,1) -1;
      stimOnsets{iStim} = stimOnsets{iStim}(stimPresent);
    end
  end
end
if ieNotDefined('runTransitions')
  runTransitions = [1 1];
  for iStim = 1:length(stimOnsets)
    if ~isempty(stimOnsets{iStim})
      runTransitions(2) = max(runTransitions(2),stimOnsets{iStim}(end));
    end
  end
end
 

%fill the stim matrix
stimMatrix = zeros(runTransitions(end,2)-runTransitions(1,1)+1,length(stimOnsets));
for iRun = 1:size(runTransitions,1)
  thisStimMatrix = zeros(runTransitions(iRun,2)-runTransitions(iRun,1)+1,length(stimOnsets));
  for iStim = 1:length(stimOnsets)
    ons = stimOnsets{iStim}(stimOnsets{iStim}>=runTransitions(iRun,1) & ...
      stimOnsets{iStim}<=runTransitions(iRun,2) ) - ...
      runTransitions(iRun,1)+1;
    % detect non-integer onset indices. This can happen when doing
    % multi-shot EPI and the experiment misses the volume acquisition
    % trigger and catches a non-volume shot trigger instead. So we floor
    % because triggers are more likely to be late than early
    roundons = floor(ons);
    if ~all(ons == roundons)
      % make sure the rounding hasn't interpolated onsets outside the run
      roundons(roundons<runTransitions(iRun,1)) = runTransitions(iRun,1);
      roundons(roundons>runTransitions(iRun,2)) = runTransitions(iRun,2);
      inter = roundons - ons;
      warning(['rounding onsets to nearest full volume ' ...
        '(interpolating min=%.2f, max=%.2f mean=%.2f volumes)'],...
        min(inter),max(inter),mean(inter));
      ons = roundons;
    end
    %only use stimOnsets{iStim}s that are within this runs volume numbers
    thisStimMatrix(ons,iStim ) = 1;
  end
  stimMatrix(runTransitions(iRun,1):runTransitions(iRun,2),:)=thisStimMatrix;
end
