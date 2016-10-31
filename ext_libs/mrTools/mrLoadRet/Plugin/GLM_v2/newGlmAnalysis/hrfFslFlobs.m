% hrfFslFlobs.m
%
%        $Id: hrfFslFlobs.m 1950 2010-12-18 10:12:48Z julien $
%      usage: [params,hrf] = hrfFslFlobs(params, sampleDuration, notUsed, defaultParams)
%         by: by julien besle
%       date: 14/06/07, 09/02/2010
%    purpose: reads a basis set from a flobs file
%
function [params, hrf] = hrfFslFlobs(params, sampleDuration, sampleDelay, defaultParams)

if ~any(nargin == [1 2 3 4])% 5])
  help hrfDoubleGamma
  return
end

fslPath = mrGetPref('fslPath');
if strcmp(mrGetPref('fslPath'),'FSL not installed')
  mrErrorDlg('(applyFslTFCE) No path was provided for FSL. Please set MR preference ''fslPath'' by running mrSetPref(''fslPath'',''yourpath'')')
end


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
  params.description = 'Flobs basis set';
end
if fieldIsNotDefined(params,'flobsBasisSetFile')
  params.flobsBasisSetFile = [fslPath(1:end-4) '/etc/default_flobs.flobs/hrfbasisfns.txt'];
  if ~ exist(params.flobsBasisSetFile,'file')
    mrWarnDlg(sprintf('(hrfFslFlobs) Could not find default FLOBS file ''%s''. Please check your fslPath preference',params.flobsBasisSetFile));
  end
end
  
paramsInfo = {...
    {'description', params.description, 'Comment describing the HRF model'},...
    {'flobsBasisSetFile', params.flobsBasisSetFile, 'callback',@testFileExists},...
    {'chooseFile', 0, 'type=pushbutton','callback',{@getBasisSetFile,fslPath},'buttonString=Choose File...'},...
    {'makeFlobs', 0, 'type=pushbutton','callback',{@launchMakeFlobs,fslPath},'buttonString=Launch FLOBS'},...
    {'displayHRF',0,'type=pushbutton','callback',@dispModelHRF,'buttonString=Display HRF','passParams=1','callbackArg',{threshold sampleDuration sampleDelay},'Display the hrf with current parameters'},...
};
      
if defaultParams
  params = mrParamsDefault(paramsInfo);
else
  params = mrParamsDialog(paramsInfo, 'Set model HRF parameters');
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

modelHrf = dlmread(params.flobsBasisSetFile);
%the sampling rate of Flobs basis sets is 20Hz
dt = 0.05;
t = dt*(0:(size(modelHrf,1)-1));

%normalise so that integral of sum = 1
modelHrf = modelHrf./sum(modelHrf(:));
    
%downsample with constant integral
dsFactor = round(sampleDuration/dt);
dsDelay = floor(rem(sampleDelay,sampleDuration)/dt)+1;
hrf = mrDownsample(modelHrf, dsFactor, dsDelay);

%remove trailing zeros
%hrf = hrf(1:end-find(flipud(max(abs(hrf),[],2))>threshold,1,'first')+1,:);

%output the max amplitude of the actual model HRF
params.maxModelHrf = sampleDuration/dt * max(modelHrf); 

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

if ~ exist(params.flobsBasisSetFile,'file')
    mrWarnDlg('(hrfFslFlobs) Please load a valid FLOBS output file');
    return
end

global modelHRFFig;

[params hrf t] = getModelHrf(params,callbackArg{1},callbackArg{2},callbackArg{3});
modelHRFFig = mlrSmartfig('hrfFlobs','reuse');clf;
plot(t,hrf,'o-','linewidth',2);
title({'Model HRF (under the assumption that points are sampled','in the middle of the frame period. This can be changed in the next menu)'});
xlabel('Time (sec)');
ylabel('Response magnitude');


function launchMakeFlobs(fslPath)

if ismac
  flobsCommand = 'Make_flobs_gui';
else %on linux (at least Ubuntu), the command is Make_flobs
  flobsCommand = 'Make_flobs';
end

try
  [s,w] = unix(sprintf('%s/%s',fslPath,flobsCommand));
  if s ~=- 0 % unix error
    disp('UNIX error message:')
    disp(w)
    disp('-------------------')
    return
  end
catch 
  disp(sprintf('(hrfFslFlobs) There was a problem running the %s unix command',flobsCommand))
  disp(sprintf('unix error code: %d; %s', s, w))
  return
end

function testFileExists(params)

if ~exist(params.flobsBasisSetFile,'file')
  mrWarnDlg([params.flobsBasisSetFile ' does not exist']);
end


function getBasisSetFile(fslPath)

[basisSetFilename pathname] = uigetfile('*.*','FLOBS Basis Set Text File');
if isnumeric(basisSetFilename)
   return
end

params.flobsBasisSetFile = [pathname basisSetFilename];

mrParamsSet(params);

