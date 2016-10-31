% getGlmAdvancedStatisticParamsGUI.m
%
%        $Id: getGlmAdvancedStatisticParamsGUI.m 2013 2011-01-21 17:43:19Z julien $
%      usage: params = getGlmAdvancedStatisticParamsGUI(params,useDefault)
%         by: julien besle, 
%       date: 19/01/2011
%    purpose: returns additional statistical test parameters for GLM analysis
%

function params = getGlmAdvancedStatisticParamsGUI(thisView,params,useDefault)

keepAsking = 1;

while keepAsking
  %Default params
  if fieldIsNotDefined(params,'covCorrectionMethod')
     params.covCorrectionMethod = 'generalizedLeastSquares';
  end
  if fieldIsNotDefined(params,'covEstimation')
     params.covEstimation = 'singleTukeyTapers';
  end
  if fieldIsNotDefined(params,'covFactorization')
     params.covFactorization = 'Cholesky';
  end
  if fieldIsNotDefined(params,'outputStatistic')
    params.outputStatistic = 0;
  end
  if fieldIsNotDefined(params, 'testOutput')
      params.testOutput = mrGetPref('statisticalTestOutput');
  end
  if fieldIsNotDefined(params, 'bootstrapIntervals')
      params.bootstrapIntervals = 0;
  end
  if fieldIsNotDefined(params, 'alphaConfidenceIntervals')
      params.alphaConfidenceIntervals = 0.05;
  end
  if fieldIsNotDefined(params, 'TFCE')
    params.TFCE = 0;
  end
  if fieldIsNotDefined(params,'parametricTests')
    params.parametricTests = 1;
  end
  if fieldIsNotDefined(params,'bootstrapTests')
    params.bootstrapTests = 0;
  end
  if fieldIsNotDefined(params,'permutationTests')
    params.permutationTests = 0;
  end
  if fieldIsNotDefined(params, 'nResamples')
    params.nResamples = 10000;
  end
  if fieldIsNotDefined(params, 'bootstrapFweAdjustment')
    params.bootstrapFweAdjustment = 0;
  end
  if fieldIsNotDefined(params, 'permutationFweAdjustment')
    params.permutationFweAdjustment = 0;
  end
  if fieldIsNotDefined(params, 'resampleFweMethod')
    params.resampleFweMethod = 'Adaptive Step-down';
  end
  if fieldIsNotDefined(params, 'fweAdjustment')
      params.fweAdjustment = 0;
  end
  if fieldIsNotDefined(params, 'fweMethod')
      params.fweMethod = 'Adaptive Step-down';
  end
  
  if fieldIsNotDefined(params, 'fdrAdjustment')
      params.fdrAdjustment = 0;
  end
  if fieldIsNotDefined(params, 'fdrAssumption')
      params.fdrAssumption = 'Independence/Positive dependence';
  end
  if fieldIsNotDefined(params, 'fdrMethod')
      params.fdrMethod = 'Adaptive Step-up';
  end
  
  if fieldIsNotDefined(params, 'trueNullsEstimationMethod')
      params.trueNullsEstimationMethod = 'Least Squares';
  end
  if fieldIsNotDefined(params, 'trueNullsEstimationThreshold')
      params.trueNullsEstimationThreshold = .05;
  end
  if fieldIsNotDefined(params, 'noDoubleBootstrap')
      params.noDoubleBootstrap = 0;
  end
  
  if strcmp(mrGetPref('fslPath'),'FSL not installed')
    params.TFCE = 0;
    tfceOptionVisible = 'enable=0';
  else
      tfceOptionVisible = 'visible=1';
  end

  covCorrectionMethodMenu = putOnTopOfList(params.covCorrectionMethod,{'varianceCorrection','preWhitening','generalizedLeastSquares'});%, 'generalizedFTest'});
  covEstimationMenu = putOnTopOfList(params.covEstimation,{'singleTukeyTapers','dampenedOscillator'});
  covFactorizationMenu = putOnTopOfList(params.covFactorization,{'Cholesky'});
  testOutputMenu = putOnTopOfList(params.testOutput,{'P value','Z value','-log10(P) value'});
  resampleFweMethodMenu = putOnTopOfList(params.resampleFweMethod,{'Single-step','Adaptive Single-step','Step-down','Adaptive Step-down'});
  fdrAssumptionMenu = putOnTopOfList(params.fdrAssumption,{'Independence/Positive dependence','None'});
  fweMethodMenu = putOnTopOfList(params.fweMethod,{'Single-step (Bonferroni)','Adaptive Single-step','Step-down (Holm)','Adaptive Step-down','Step-up (Hochberg)','Hommel'});
  fdrMethodMenu = putOnTopOfList(params.fdrMethod,{'Step-up','Adaptive Step-up','Two-stage Step-up','Multiple-stage Step-down','Multiple-stage Step-up'});
  trueNullsEstimationMethodMenu = putOnTopOfList(params.trueNullsEstimationMethod,{'Lowest Slope','Raw P-values Cut-off','Least Squares'});
  
  paramsInfo = {...
    {'covCorrection',params.covCorrection,'type=checkbox','(EXPERIMENTAL) Correction for temporally-correlated noise. Correcting for noise correlation is important for single-subject level statistics but significantly increases analysis time. Uncorrected correlated noise biases statistical significance of contrasts/F-tests but should not affect parameter estimates (contrasts values).'},...
    {'covCorrectionMethod',covCorrectionMethodMenu,'type=popupmenu','contingent=covCorrection','Type of noise correlation correction. Generalized least square is the preferred method, but pre-whitening gives identical results and is more efficient in conjunction with bootstrap/permutation statistics. Variance correction is less conservative but faster (see Wiki for the different algorithms and computing time comparison).'},...
    {'covEstimation',covEstimationMenu,'visible=0','type=popupmenu','contingent=covCorrection','Type of Estimation of the noise covariance matrix'},...
    {'covFactorization',covFactorizationMenu,'visible=0','type=popupmenu','contingent=covCorrection','Type of factorization of the covariance matrix for the computation of the pre-whitening filter'},...
    {'outputStatistic', params.outputStatistic,'type=checkbox', 'Outputs the statistic (T/F) as an overlay in addition to the Z-transform/P value'},...
    {'testOutput', testOutputMenu,'type=popupmenu', 'Type of statistics for output overlay.  P: outputs the probability value associated with the T/F statistic. Z: outputs standard normal values associated with probability p. Default output can be set using mrSetPref. P-values less than 1e-16 will be replaced by 0, corresponding Z values by +/-8.209536145151493, and corresponding -10log(P) values by 15.653559774527022'},...
    {'parametricTests', params.parametricTests,'type=checkbox', 'Performs parametric T-tests on contrasts and parametric F-tests. Disable this if you only want bootstrap or permutation-based T/F-tests.'},...
    {'fweAdjustment', params.fweAdjustment,'type=checkbox', 'Whether to perform adjustments for Bonferroni-type familywise error rate (FWER) control'},...
    {'fweMethod', fweMethodMenu,'contingent=fweAdjustment','type=popupmenu', 'Bonferroni-type FWER control method. Hommel method is the most powerful of non-adaptive methods. Adaptive methods estimate the number of true Null hypotheses from the raw p-values, making the procedure more powerful, but are only implemented for single-step and step-up methods.'},...
    {'fdrAdjustment', params.fdrAdjustment,'type=checkbox', 'Whether to perform adjustments for false discovery rate (FDR) control'},...
    {'fdrAssumption', fdrAssumptionMenu,'contingent=fdrAdjustment','type=popupmenu', 'Distributional assumption for the FDR adjustment methods. Most FDR methods assume that tests are independent or positively correlated, although some have a correcting factor in case no such assumption is made'},...
    {'fdrMethod', fdrMethodMenu,'contingent=fdrAdjustment','type=popupmenu', 'Type of FDR control adjustment method. The multistep adaptive method is supposedly the most powerful'},...
    {'trueNullsEstimationMethod', trueNullsEstimationMethodMenu,'type=popupmenu', 'Method to estimate the number of true null outcomes for adaptive adjustment methods'},...
    {'trueNullsEstimationThreshold', params.trueNullsEstimationThreshold, 'Raw p-values cutoff for true nulls estimation'},...
    {'bootstrapTests', params.bootstrapTests,'type=checkbox', '(EXPERIMENTAL) Performs non-parametric residual bootstrap tests on T/F values. Bootstrapping consists in resampling the residuals with replacement after OLS/GLS fit, using these bootstrapped residuals as the new time-series for OLS/GLS fitting and estimating the null-hypothesis distributions for nResamples resamplings.'},...
    {'permutationTests', params.permutationTests,'type=checkbox', '(EXPERIMENTAL) Performs non-parametric permutation tests on contrasts/F values. Permutations constist in shuffling stimulus event labels, re-fitting the GLM using OLS/GLS and estimating the null-hypothesis distribution for each statistic for nResamples permutations.'},...
    {'nResamples', params.nResamples, 'minmax=[10 inf]', 'Number of iterations for bootstrap/permutation tests.'},...
    {'bootstrapIntervals', params.bootstrapIntervals,'type=checkbox', '(EXPERIMENTAL) Whether to compute bootstrap confidence intervals for contrasts and parameter estimates'},...
    {'alphaConfidenceIntervals', params.alphaConfidenceIntervals,'contingent=bootstrapIntervals', 'minmax=[0 1]', 'Confidence Intervals will be computed as quantiles [alpha/2 1-alpha/2] of the bootstrap-estimated null distribution'},...
    {'bootstrapFweAdjustment', params.bootstrapFweAdjustment,'type=checkbox', '(EXPERIMENTAL) Adjusts P-values for familywise error rate control by residuals bootstrapping.'},...
    {'permutationFweAdjustment', params.permutationFweAdjustment,'type=checkbox', '(EXPERIMENTAL) Adjusts P-values for familywise error rate control by permutations.'},...
    {'resampleFweMethod', resampleFweMethodMenu,'type=popupmenu', 'Bootstrap/permutation-based family-wise error control method. Adaptive methods estimate the number of true Null hypotheses form the raw p-values. Adaptive Step-down is the most powerful'},...
    {'noDoubleBootstrap', params.noDoubleBootstrap,'type=checkbox', 'Prevents boostrap-based FWER control adjustment on bootstrap P-value to spare memory'},...
    {'TFCE', params.TFCE,tfceOptionVisible,'type=checkbox', 'Performs Threshold Free Cluster Enhancement on T/F maps using fslmaths. This option is only enabled if a path is specified for FSL by running mrSetPref(''fslPath'',''yourpath''). In addition it can only be used in conjonction with bootstrap/permutation tests or resample-based FWER control adjustment.'},...
       };

  if ~params.showAdvancedStatisticMenu || useDefault
    tempParams = mrParamsDefault(paramsInfo);
  else
    tempParams = mrParamsDialog(paramsInfo,'Advanced Statistics Menu');
  end

  % user hit cancel
  if isempty(tempParams)
    params = tempParams;
    return;
  end
  
  params = mrParamsCopyFields(tempParams,params);
  %this is because of the incoherent behaviour of mrParamsGet that empties disabled params fields
  if isempty(params.TFCE)
    params.TFCE = 0;
  end
  
  %check consistency of parameters
  if (params.computeTtests || params.numberFtests) && params.permutationTests && ismember(params.scanParams{params.scanNum(1)}.stimDurationMode,{'From file','Block'})
    mrWarnDlg('(getTestParamsGUI) Permutation tests are not compatible with stimulus duration from log file or block mode','Yes');
  elseif params.TFCE && params.parametricTests && (params.bootstrapFweAdjustment || params.permutationFweAdjustment) && ismember(params.resampleFweMethod,{'Step-down','Adaptive Step-down'})
    mrWarnDlg('(getTestParamsGUI) Step-down resample-based FWE adjustment is not implemented for TFCE-transformed data','Yes');
  elseif params.TFCE && params.parametricTests && (params.bootstrapFweAdjustment || params.permutationFweAdjustment) && ismember(params.resampleFweMethod,{'Adaptive Single-step'})
    mrWarnDlg('(getTestParamsGUI) Adaptive resample-based FWE adjustment is not implemented for TFCE-transformed data','Yes');
  elseif params.fdrAdjustment && strcmp(params.fdrAssumption,'None') && ismember(params.fdrMethod,{'Adaptive Step-down','Multiple-stage Adaptive Step-up'})
    mrWarnDlg('(getTestParamsGUI) Multi-stage and Step-down adaptive FDR adjustments require the assumption that voxel are independent or positively correlated','Yes');
  elseif params.bootstrapTests  && ~ismember(params.analysisVolume,{'Loaded ROI(s)' 'Visible ROI(s)'})
    mrWarnDlg('(getTestParamsGUI) Bootstrap tests are currently only allowed for ROI(s) analyses','Yes');
  elseif params.bootstrapFweAdjustment && params.bootstrapTests  && ~ismember(params.analysisVolume,{'Loaded ROI(s)' 'Visible ROI(s)'})
    mrWarnDlg('(getTestParamsGUI) Resampled-based FWE adjustments are currently only allowed for ROI(s) analyses','Yes');
  elseif params.permutationTests && ~all( (sum(logical(allContrasts),2)==2 & sum(allContrasts,2)==0) | ~any(allContrasts,2))
    mrWarnDlg('(getTestParamsGUI) Randomization tests can only be run if all contrasts of all T/F tests are 1 to 1 comparisons','Yes');
  elseif params.TFCE && ~(params.bootstrapTests || params.bootstrapFweAdjustment || params.permutationTests || params.permutationFweAdjustment)
    mrWarnDlg('(getTestParamsGUI) TFCE requires resampling tests (bootstrap, permutation or resample-based FWE control)','Yes');
  elseif params.permutationTests && params.permutationFweAdjustment && ~params.parametricTests
    mrWarnDlg('(getTestParamsGUI) permutation-based adjustment for permutation tests is not implemented','Yes');
  elseif params.bootstrapFweAdjustment && ~(params.bootstrapTests ||params.parametricTests)
    mrWarnDlg('(getTestParamsGUI) bootstrap-based FWE adjustments must be performed on bootsrap or parametric tests','Yes');
  elseif params.permutationFweAdjustment && ~(params.permutationTests ||params.parametricTests)
    mrWarnDlg('(getTestParamsGUI) permutation-based FWE adjustments must be performed on permutation or parametric tests','Yes');
  elseif (params.fweAdjustment || params.fdrAdjustment) && ~(params.permutationTests ||params.parametricTests||params.bootstrapTests)
    mrWarnDlg('(getTestParamsGUI) FWE/FDR adjustments must be performed on bootstrap/permutation/parametric tests','Yes');
  else
    keepAsking = 0;
  end
  if keepAsking && useDefault
    params = [];
    return;
  end

end


