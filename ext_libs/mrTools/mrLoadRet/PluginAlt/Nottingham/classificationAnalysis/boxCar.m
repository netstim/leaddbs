% boxCar.m
%
%        $Id: boxCar.m 2081 2011-03-04 10:46:48Z julien $
%      usage: [params,hrf] = boxCar(params, sampleDuration, stimDuration )
%         by: julien besle
%       date: 14/06/07, 09/02/2010
%    purpose: returns a canonical hrf modeled as a boxcar function
%
function [params, hrf] = boxCar(params, sampleDuration, notUsed, defaultParams )

threshold = 1e-3; %threshold for removing trailing zeros at the end of the model

if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('sampleDuration'),sampleDuration = 1;end

if ieNotDefined('params')
  params = struct;
end
if ~isfield(params,'description')
  params.description = 'Box Car Function';
end
if ~isfield(params,'hdLag')
  params.hdLag = 2;
end
  
paramsInfo = {...
    {'description', params.description, 'Comment describing the hdr model'},...
    {'hdLag', params.hdLag, 'Estimate of the HD lag, shifts the boxcar by this number of TRs'},...
};
      
if defaultParams
  params = mrParamsDefault(paramsInfo);
else
  params = mrParamsDialog(paramsInfo, 'Set model HRF parameters');
end

if nargout==1
   return
end


hrf = [zeros(1,params.hdLag),1]';




% 
% 
% 
% 
% 
% tmax = max(params.y*3, 20); %min length of the hrf model in seconds
% 
% if isfield(params, 'tmax')
%     tmax = params.tmax;
% end
% 
% shift = 0;
% if isfield(params, 'shift')
%     shift = params.shift;
% end
% 
% dt = 0.05;
% 
% t = 0:dt:tmax;
% warning('off', 'MATLAB:log:logOfZero');
% modelHrf = gampdf(t, params.x, 1) - gampdf(t, params.y, 1)/params.z;
% warning('on', 'MATLAB:log:logOfZero');
% 
% if shift<0
%   modelHrf = [zeros(1, ceil(-shift/dt)), modelHrf];
% elseif shift>0
%   modelHrf = modelHrf( ceil(shift/dt):end );
% end
% 
% 
% if params.includeDerivative
%   % take the derivative
%   modelHrfDerivative = [diff(modelHrf), 0];
%   % orthogonalize
%   modelHrfDerivative = modelHrfDerivative - modelHrf*(modelHrfDerivative/modelHrf);
%   % remove mean
%   modelHrfDerivative = modelHrfDerivative - mean(modelHrfDerivative);
%   % normalize so that it's norm equals the Hrf norm
%   modelHrfDerivative = modelHrfDerivative / norm(modelHrfDerivative)*norm(modelHrf);
%   %concatenate
%   modelHrf = [modelHrf; modelHrfDerivative];
% end
% 
% %remove trailing zeros
% modelHrf = modelHrf(1:end-find(flipud(max(abs(modelHrf),[],2))>threshold,1,'first')+1,:);
% %normalise so that integral of sum = 1
% modelHrf = modelHrf./sum(modelHrf(:));
%     
% %downsample with constant integral
% hrf = downsample(modelHrf', round(sampleDuration/dt));
% 
% 
% params.maxModelHrf = sampleDuration/dt * max(modelHrf'); %output the max amplitude of the actual model HRF
% 
