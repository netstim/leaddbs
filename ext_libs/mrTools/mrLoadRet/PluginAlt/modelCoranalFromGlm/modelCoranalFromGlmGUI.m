% modelCoranalFromGlmGUI
%
%        $Id$
%      usage: params = modelCoranalFromGlmGUI('thisView',thisView,'groupName',groupName,'params',params)
%         by: julien besle
%       date: 2010-10-24
%     inputs: 
%    outputs: 
% 
%    purpose: parameters for modelCoranalFromGlm.m
%        e.g:
%

function params = modelCoranalFromGlmGUI(varargin)

% get the arguments
eval(evalargs(varargin));
if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('thisView'),thisView = newView;end
if ieNotDefined('groupName'), groupName = viewGet(thisView,'groupName');end;

if ieNotDefined('params')
  params = struct;
  if isempty(thisView.analyses) || ...
      ~ismember(thisView.analyses{thisView.curAnalysis}.type,{'glmAnalStats','glmAnal','erAnal','deconvAnal'})
    mrWarnDlg('You must load an event-related/GLM Analysis');
    params = [];
    return;
  end
  
  analysisParams = convertOldGlmParams(viewGet(thisView,'analysisParams'));
  params.originalAnalysis = analysisParams.saveName;
  params.originalGroup = analysisParams.groupName;

  d = viewGet(thisView,'d');
  if fieldIsNotDefined(params,'newGroupName')
    params.newGroupName = ['Average modeled from ' groupName ' GLM'];
  end
  if fieldIsNotDefined(params,'saveName')
    params.saveName = 'corAnalFromGLM';
  end

  if fieldIsNotDefined(params,'ncycles')
    params.ncycles = 12;
  end
  if fieldIsNotDefined(params,'stimDurationS')
    params.stimDurationS = 2*d.tr;
  end
  if fieldIsNotDefined(params,'cycleLengthTR')
    params.cycleLengthTR = 10;
  end
  if fieldIsNotDefined(params,'dummyStimsTR')
    params.dummyStimsTR = 6;
  end
  if fieldIsNotDefined(params,'numberStims')
    params.numberStims = d.nhdr;
  end
  if fieldIsNotDefined(params,'nScans')
    %this is the maximum even number of scans we can have
    params.nScans = 2*floor(size(d.scm,1)/(params.ncycles*params.cycleLengthTR+params.dummyStimsTR)/2);
  end
  if fieldIsNotDefined(params,'reverseScans')
    params.reverseScans = mat2str(2:2:params.nScans);
  end
  if fieldIsNotDefined(params,'delayTR')
    params.delayTR = -2;
  end
  if fieldIsNotDefined(params,'subsetBox')
    params.subsetBox = mat2str([1 size(d.ehdr,1);1 size(d.ehdr,2);1 size(d.ehdr,3)]);
  end
  if fieldIsNotDefined(params,'TRCorrection')
    params.TRCorrection = 'Yes';
  end
  if fieldIsNotDefined(params,'noiseMode')
    params.noiseMode = 'Residuals';
  end
  if fieldIsNotDefined(params,'detrend')
    params.detrend = 'Highpass';
  end
  if fieldIsNotDefined(params,'spatialnorm')
    params.spatialnorm = 'None';
  end
  if fieldIsNotDefined(params,'trigonometricFunction')
    params.trigonometricFunction = 'Cosine';
  end
elseif isfield(params,'scanParams')
  params = params.scanParams{viewGet(thisView,'curScan')};
  %if params is not empty, we assume that they are all present and correct
  if fieldIsNotDefined(params,'trigonometricFunction') %except for this one, which has been added later
    params.spatialnorm = 'Cosine';
  end
end

groupNames = putOnTopOfList(params.originalGroup,viewGet(thisView,'groupNames'));
TRCorrectionMenu = putOnTopOfList(params.TRCorrection,{'Yes','No'});
noiseModeMenu = putOnTopOfList(params.noiseMode,{'Residuals','No noise'});
detrendMenu = putOnTopOfList(params.detrend,{'Highpass','None','Linear','Quadratic'});
spatialnormMenu = putOnTopOfList(params.spatialnorm,{'Divide by mean','None'});
trigonometricFunctionMenu = putOnTopOfList(params.trigonometricFunction,{'Sine','Cosine'});

askForParams = 1;
while askForParams
  paramsInfo = {...
      {'originalGroup',groupNames,'type=popupmenu','Name of the Group from which the GLM anlaysis is taken'},...
      {'originalAnalysis',params.originalAnalysis,'Name of the GLM analysis from which the estimates are taken'},...
      {'newGroupName',params.newGroupName,'Name of the group in which the created scans/analysis will be install'},...
      {'saveName',params.saveName,'Name of the analysis to save. If identical to an existing analysis a scan will be added and overlays will be added to the existed analysis'},...
      {'numberStims',params.numberStims,'number of different stimulations, each to be modelled as a linear contrast of GLM estimates, defined in a following dialog box'},...
      {'nScans',params.nScans,'Number of scans to model'},...
      {'reverseScans', params.reverseScans, 'Index of scans with reverse sequence'},...
      {'ncycles', params.ncycles, 'number of stimulation cycles per scan'},...
      {'cycleLengthTR', params.cycleLengthTR, 'Length of a stimulation cycle in TRs'},...
      {'stimDurationS', params.stimDurationS, 'Duration of each stimulation'},...
      {'dummyStimsTR', params.dummyStimsTR,  'number of TR separating the start of the stimulation and the start of the scan acquisition'},...
      {'delayTR', params.delayTR, 'By how many TRs to delay the timeseries before avaraging'},...
      {'TRCorrection', TRCorrectionMenu, 'type=popupmenu', 'Whether to delay reversed stimulation time-series by an extra TR to correct for the fact that all acquisitions are separated by one TR'},...
      {'noiseMode', noiseModeMenu, 'type=popupmenu', 'Residuals uses the residuals of the GLM analysis as an estimate for the noise. If no noise is required, only two scans are modelled because the result is indpendent of the number of scans.'},...
      {'detrend', detrendMenu, 'type=popupmenu', 'Whether to detrend the raw time-series'},...
      {'spatialnorm', params.spatialnorm, 'editable=0', 'Spatial normalization is disabled because time-series mean might be negative'},...
      {'trigonometricFunction', trigonometricFunctionMenu, 'type=popupmenu', 'Whether to take the phase of the sine or cosine functions'},...
      {'subsetBox', params.subsetBox, 'subset of voxels, of the form [X1 X2;Y1 Y2;Z1 Z2] (Zs are optional)'},...
    };

  if defaultParams
    params = mrParamsDefault(paramsInfo);
  else
    params = mrParamsDialog(paramsInfo, 'Correlation from GLM parameters');
  end
  % Abort if params empty
  if ieNotDefined('params'),return,end

  groupNames = putOnTopOfList(params.originalGroup,groupNames);
  TRCorrectionMenu = putOnTopOfList(params.TRCorrection,TRCorrectionMenu);
  noiseModeMenu = putOnTopOfList(params.noiseMode,noiseModeMenu);
  detrendMenu = putOnTopOfList(params.detrend,detrendMenu);
  spatialnormMenu = putOnTopOfList(params.spatialnorm,spatialnormMenu);
  trigonometricFunctionMenu = putOnTopOfList(params.trigonometricFunction,trigonometricFunctionMenu);

  if ~ismember(params.originalAnalysis,viewGet(thisView,'analysisNames',viewGet(thisView,'groupNum',groupNames{1})))
    mrWarnDlg(['(modelCoranalFromGlmGUI) There is no analysis ' params.originalAnalysis ' in group ' groupNames{1}] );
  elseif max(eval(params.reverseScans))>params.nScans
    mrWarnDlg('(modelCoranalFromGlmGUI) Max(reverseScan) must be less than number of scans');
  elseif params.nScans/length(eval(params.reverseScans))~=2
    mrWarnDlg('(modelCoranalFromGlmGUI) The number of reverse scans must equal half the total number of scans');
  elseif ~ieNotDefined('d') &&...
      strcmp(params.noiseMode,'Residuals') && params.nScans*(params.ncycles*params.cycleLengthTR+params.dummyStimsTR) >size(d.scm,1)
      mrWarnDlg('(modelCoranalFromGlmGUI) If residuals are used as noise, the total number of simulated TRs cannot be greater than the number of TRs in the GLm analysis');
  elseif ~ieNotDefined('d') && ...
      rem(params.stimDurationS,d.tr) || rem(params.cycleLengthTR,params.numberStims)
    mrWarnDlg('(modelCoranalFromGlmGUI) Function not implemented for stimulus durations that are not multiples of TR or cycle lengths that are not multiples of the number of stimuli/contrasts');
  else
    while askForParams
      
      tempParams = getEstimateContrastsGUI(params,analysisParams,defaultParams);
      
      if isempty(tempParams)
        if size(tempParams,2)
          params = tempParams;
          return;
        else
          break;
        end
      else
        askForParams=0;
        params = tempParams;
      end
    end
  end
end
               

function params = getEstimateContrastsGUI(params,analysisParams,useDefault)

if fieldIsNotDefined(params, 'contrasts') || ...
  ~isequal(size(params.contrasts,1),params.numberStims) || ~isequal(size(params.contrasts,2),analysisParams.numberEVs)
    contrasts = eye(params.numberStims,analysisParams.numberEVs);
else
  if ischar(params.contrasts)
    contrasts = eval(params.contrasts);
  else
    contrasts = params.contrasts;
  end
end

if fieldIsNotDefined(analysisParams,'EVnames')
  d = viewGet(thisView,'d');
  EVnames = repmat({' '},1,length(d.stimNames));
  if length(d.stimNames) < params.numberStims
    EVnames = d.stimNames(1:analysisParams.numberEVs);
  else
    EVnames(1:params.numberStims) = d.stimNames;
  end
else
  EVnames = analysisParams.EVnames;
end

%shift the order of contrasts to account for dummies
contrasts = circshift(contrasts,round(params.dummyStimsTR/params.cycleLengthTR*size(contrasts,1)));

paramsInfo = cell(1,params.numberStims+1);
paramsInfo{1}= {'EVnames', EVnames, 'type=stringarray','editable=0','Name of the estimated EVs in the GLM analysis'};
for iStim = 1:params.numberStims
  paramsInfo{iStim+1} = {['Stimulus_' num2str(iStim)], contrasts(iStim,:),...
    'incdec=[-1 1]','incdecType=plusMinus','minmax=[0 inf]',...
    ['Combination of GLM estimates defining stimulus ' num2str(iStim)]};
end

%%%%%%%%%%%%%%%%%%%%%%%
% now we have all the dialog information, ask the user to set parameters
if useDefault
   tempParams = mrParamsDefault(paramsInfo);
else
   tempParams = mrParamsDialog(paramsInfo,'Set Stimulus/EV matrix');
end

% user hit cancel
if isempty(tempParams)
  params=tempParams;
  return;
end

% form contrats matrix from fields
contrasts = zeros(params.numberStims,length(EVnames));
for iStim = 1:params.numberStims
  contrasts(iStim,:) = tempParams.(['Stimulus_' num2str(iStim)]);
  tempParams = mrParamsRemoveField(tempParams,['Stimulus_' num2str(iStim)]);
end
params = mrParamsCopyFields(...
  mrParamsDefault({{'contrasts',contrasts,'Matrix forming stimuli from combinations of GLM estimates'}}),params);

%copy the EV names
params = mrParamsCopyFields(tempParams,params);
