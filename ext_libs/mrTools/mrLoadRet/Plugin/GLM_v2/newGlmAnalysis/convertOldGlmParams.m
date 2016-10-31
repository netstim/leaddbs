% convertOldGlmParams.m
%
%      usage: params = convertOldGlmParams(params)
%         by: julien besle 
%       date: 04/12/2010
%    purpose: re-organizes old GLM params for new GLM code
%              $Id$

function params = convertOldGlmParams(params)

if isfield(params,'testParams')
  params = copyFields(params.testParams,params);
  params = rmfield(params,'testParams');
end

if isfield(params,'contrast') 
  if ~isfield(params,'contrasts')
    params.contrasts = params.contrast;
  end
  params = rmfield(params,'contrast');
end

if isfield(params,'f_tests')
  if ischar(params.f_tests)
    params.fTests = eval(params.f_tests);
  else
    params.fTests = params.f_tests;
  end
  params.numberFtests = size(params.fTests,1);
  params = rmfield(params,'f_tests');
end
if isfield(params,'params')  
  if isfield(params,'fTests')
    for iFtest = 1:size(params.fTests,1)
      params.restrictions{iFtest} = diag(params.fTests(iFtest,:));
    end
    params = rmfield(params,'fTests');
  end
end

if isfield(params,'hrfModel') && strcmp(params.hrfModel,'hrfDiffGamma')
  params.hrfModel='hrfDoubleGamma';
end
if isfield(params,'inplaceConcat') || isfield(params,'applyFiltering')
  params.hrfModel='hrfDeconvolution';
end
  
if isfield(params,'hrfParams')  
  if isfield(params.hrfParams,'incDeriv')
    params.hrfParams.includeDerivative = params.hrfParams.incDeriv;
    params.hrfParams = rmfield(params.hrfParams,'incDeriv');
  end
end

%scan-specific parameters
if isfield(params,'scanParams')
  for iScan = 1:length(params.scanParams)
    if ~isempty(params.scanParams{iScan})
      if isfield(params.scanParams{iScan},'trSupersampling')
        params.scanParams{iScan}.designSupersampling = params.scanParams{iScan}.trSupersampling;
        params.scanParams{iScan} = rmfield(params.scanParams{iScan},'trSupersampling');
      end
      if ~isfield(params.scanParams{iScan},'supersamplingMode')
        params.scanParams{iScan}.supersamplingMode='Automatic';
      end
      if isfield(params.scanParams{iScan},'forceStimOnSampleOnset') && params.scanParams{iScan}.forceStimOnSampleOnset
        params.scanParams{iScan}.supersamplingMode = 'Set value';
        params.scanParams{iScan} = rmfield(params.scanParams{iScan},'forceStimOnSampleOnset');
      end
      if isfield(params.scanParams{iScan},'stimDuration') && ~isfield(params.scanParams{iScan},'stimDurationMode')
        if strcmp(params.scanParams{iScan}.stimDuration,'fromFile')
          params.scanParams{iScan}.stimDurationMode = 'From file';
          params.scanParams{iScan}.stimDuration = [];
        else
          params.scanParams{iScan}.stimDurationMode = 'Set value';
        end
      end
    end
  end
end

if isfield(params,'correctionType')
  params.covCorrectionMethod = params.correctionType;
  params = rmfield(params,'correctionType');
end
if isfield(params,'n_rand')
  params.nResamples = params.n_rand;
  params = rmfield(params,'n_rand');
end
if isfield(params,'outputZStatistic')
  params = rmfield(params,'outputZStatistic');
end
if isfield(params,'outputPValue')
  params = rmfield(params,'outputPValue');
end
if isfield(params,'parametricTestOutput')
  if strcmp(params.parametricTestOutput,'T/F value')
    params.outputParametricStatistic =1;
  else
    params.testOutput = params.parametricTestOutput;
  end
  params = rmfield(params,'parametricTestOutput');
end
if isfield(params,'randomizationTestOutput')
  params = rmfield(params,'randomizationTestOutput');
end
if isfield(params,'bootstrapTestOutput')
  params = rmfield(params,'bootstrapTestOutput');
end
if isfield(params, 'maskContrastOverlay')
  if params.maskContrastOverlay
    params.alphaContrastOverlay = 'FDR';
  else
    params.alphaContrastOverlay = 'None';
  end
  params = rmfield(params,'maskContrastOverlay');
end
