% hrfBoxcar.m
%
%        $Id: hrfBoxcar.m 2332 2011-09-21 14:42:18Z julien $
%      usage: [params,hrf] = hrfBoxcar(params, sampleDuration, notUsed, defaultParams)
%         by: julien besle
%       date: 13/04/2010
%    purpose: returns a canonical hrf modeled as a boxcar function
%
function [params,hrf] = hrfBoxcar(params, sampleDuration, notUsed, defaultParams)

if ~any(nargin == [1 2 3 4])% 5])
  help hrfBoxcar
  return
end
if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('sampleDuration'),sampleDuration = 1;end

%estimationSampling = varargin{1};

if ieNotDefined('defaultParams'),defaultParams = 0;end

if ieNotDefined('params')
  params = struct;
end
if fieldIsNotDefined(params,'description')
  params.description = 'Boxcar';
end
if fieldIsNotDefined(params,'delayS')
  params.delayS = 0;
end
if fieldIsNotDefined(params,'durationS')
  params.durationS = sampleDuration;
end
   
paramsInfo = {...
    {'description', params.description, 'comment describing the hdr model'},...
    {'delayS',params.delayS,'Delay before start of boxcar, in seconds'},...
    {'durationS',params.durationS,'Duration of boxcar, in seconds'},...
};
      
if defaultParams
  params = mrParamsDefault(paramsInfo);
else
  params = mrParamsDialog(paramsInfo, 'Set Boxcar parameters');
end

if nargout==1
   return
end

totalDurationS = params.delayS+params.durationS;
totalDuration = round(totalDurationS/sampleDuration);  %total duration in samples
delay = round(params.delayS/sampleDuration);  %duration of delay in samples
duration = totalDuration - delay;
hrf = [zeros(delay,1);ones(duration,1)];


%hrf = hrf(:,1:round(estimationSampling/sampleDuration):end);
