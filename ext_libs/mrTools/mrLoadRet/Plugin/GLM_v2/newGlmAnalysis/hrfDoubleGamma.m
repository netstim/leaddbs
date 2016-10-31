% hrfDoubleGamma.m
%
%        $Id$
%      usage: [params,hrf] = hrfDoubleGamma(params, sampleDuration, notUsed, defaultParams)
%         by: farshad moradi, modified by julien besle
%       date: 14/06/07, 09/02/2010
%    purpose: returns a canonical hrf that's a difference of two gamma distribution function
%
function [params, hrf] = hrfDoubleGamma(params, sampleDuration, sampleDelay, defaultParams)

threshold = 1e-3; %threshold for removing trailing zeros at the end of the model

if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('sampleDuration'),sampleDuration = 1;end
if ieNotDefined('sampleDelay')
  sampleDelay=sampleDuration/2;
end

if ieNotDefined('params')
  params = struct;
end
if fieldIsNotDefined(params,'description')
  params.description = 'Double Gamma Function';
end
if fieldIsNotDefined(params,'x')
  params.x = 6;
end
if fieldIsNotDefined(params,'y')
  params.y = 16;
end
if fieldIsNotDefined(params,'z')
  params.z = 6;
end
if fieldIsNotDefined(params,'includeDerivative')
  params.includeDerivative = 0;
end

% figure handle if one is up
global modelHRFFig;
modelHRFFig = [];

paramsInfo = {...
    {'description', params.description, 'Comment describing the hdr model'},...
    {'x', params.x, 'Shape parameter of the positive Gamma distribution function; hrf = gampdf(...,x,1)-gampdf(...,y,1)/z'},...
    {'y', params.y, 'Shape parameter of the negative Gamma distribution function; hrf = gampdf(...,x,1)-gampdf(...,y,1)/z'},...
    {'z', params.z, 'Scaling factor between the positive and negative gamma componenets; hrf = gampdf(...,x,1)-gampdf(...,y,1)/z'},...
    {'includeDerivative',params.includeDerivative,'type=checkbox','Includes derivative of the hrf in the model'},...
    {'displayHRF',0,'type=pushbutton','callback',@dispModelHRF,'buttonString=Display HRF','passParams=1','callbackArg',{threshold sampleDuration sampleDelay},'Display the hrf with current parameters'},...
};
      
if defaultParams
  params = mrParamsDefault(paramsInfo);
else
  params = mrParamsDialog(paramsInfo, 'Set model HRF parameters');
end

% close figure that displays HRF if it is up
if ~isempty(modelHRFFig) &  ishghandle(modelHRFFig)
  close(modelHRFFig);
end

if nargout==1
   return
end

% get the model HRF
if ~isempty(params)
  [params hrf] = getModelHrf(params,threshold,sampleDuration,sampleDelay);
end


%%%%%%%%%%%%%%%%%%%%%
%    getModelHrf    %
%%%%%%%%%%%%%%%%%%%%%
function [params hrf t] = getModelHrf(params,threshold,sampleDuration,sampleDelay)

tmax = max(params.y*3, 20); %min length of the hrf model in seconds

if isfield(params, 'tmax')
    tmax = params.tmax;
end

shift = 0;
if isfield(params, 'shift')
    shift = params.shift;
end

dt = 0.05;

t = 0:dt:tmax;
warning('off', 'MATLAB:log:logOfZero');
modelHrf = gampdf(t, params.x, 1) - gampdf(t, params.y, 1)/params.z;
warning('on', 'MATLAB:log:logOfZero');

if shift<0
  modelHrf = [zeros(1, ceil(-shift/dt)), modelHrf];
elseif shift>0
  modelHrf = modelHrf( ceil(shift/dt):end );
end


if params.includeDerivative
  % take the derivative
  modelHrfDerivative = [diff(modelHrf), 0];
  % orthogonalize
  modelHrfDerivative = modelHrfDerivative - modelHrf*(modelHrfDerivative/modelHrf);
  % remove mean
  modelHrfDerivative = modelHrfDerivative - mean(modelHrfDerivative);
  % normalize so that its norm equals the Hrf norm
  modelHrfDerivative = modelHrfDerivative / norm(modelHrfDerivative)*norm(modelHrf);
  %concatenate
  modelHrf = [modelHrf; modelHrfDerivative];
end

%normalise so that integral of sum = 1
modelHrf = modelHrf./sum(modelHrf(:));
    
%downsample with constant integral
dsFactor = round(sampleDuration/dt);
dsDelay = floor(rem(sampleDelay,sampleDuration)/dt)+1;
hrf = mrDownsample(modelHrf', dsFactor, dsDelay);

%remove trailing zeros
%hrf = hrf(1:end-find(flipud(max(abs(hrf),[],2))>threshold,1,'first')+1,:);

%output the max amplitude of the actual model HRF
params.maxModelHrf = sampleDuration/dt * max(modelHrf'); 

% return actual time
t = t(dsDelay:dsFactor:length(t));
% make sure t is same length as hrf
if length(t) < length(hrf)
  t = [t nan(1,length(hrf)-length(t))];
elseif length(t) > length(hrf)
  t = t(1:length(hrf));
end

%%%%%%%%%%%%%%%%%%%%%%
%    dispModelHRF    %
%%%%%%%%%%%%%%%%%%%%%%
function dispModelHRF(callbackArg,params)

global modelHRFFig;

[params hrf t] = getModelHrf(params,callbackArg{1},callbackArg{2},callbackArg{3});
modelHRFFig = mlrSmartfig('hrfDoubleGamma','reuse');clf;
plot(t,hrf,'o-','linewidth',2);
title({'Model HRF (under the assumption that points are sampled','in the middle of the frame period. This can be changed in the next menu)'});
xlabel('Time (sec)');
ylabel('Response magnitude');
