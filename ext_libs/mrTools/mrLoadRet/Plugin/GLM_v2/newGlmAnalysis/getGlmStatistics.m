% getGlmStatistics.m
%
%      usage: [d, out] = getGlmStatistics(d, params, verbose, precision, actualData, computeTtests)
%         by: Julien Besle
%       date: 18/01/10
%        $Id$
%    purpose: Fits a GLM model to timeseries at each voxel of a volume
%               and computes mass-univariate T and F statistics (parametric and bootstrap)
%               with optional correction for temporal noise correlation and bootstrap-based FWE adjustment
%
%             d needs to have a stimulus convolution matrix
%    
%             references for noise correlation correction
%             Generalized Least Squares:
%                - Wicker & Fonlupt (2003) NeuroImage, 18, p589
%                - Burock & Dale (2000) Human Brain Mapping, 11, p249
%             Pre-whitening and variance correction:
%                - Woolrich et al. (2001) NeuroImage, 14(6), p1370
%             Generalized F-tests
%                - Kruggel et al. (2002) Medical image Analysis, 6, p65
%             reference for bootstrap and FWE adjustment testing
%                - Westfall, P.H., and S.S. Young. Resampling-based multiple testing. Wiley-Interscience, 1993

function [d, out] = getGlmStatistics(d, params, verbose, precision, actualData)

%DEBUG
% lastwarn('','');

approximate_value = 0;
out = struct;
%optimalXLength = [60 80 100 120 140 160 180 200];
%optimalXLength = round(120000/nFrames);  %empirical value that optimizes the speed of bootstrapping, probably depends on available memory
optimalXLength = 60;   %optimizes the OLS residual and MSS computing time


%--------------------------------------------------- DEFAULT VALUES ------------------------------------%
if ieNotDefined('verbose'),verbose = 1;end
if ieNotDefined('precision'),precision = mrGetPref('defaultPrecision');end
if ieNotDefined('actualData'),actualData = 1;end
%if ieNotDefined('outputContrast'),outputContrast = 1;end

if ieNotDefined('params')
  params = struct;
end

%default parameters

% check that contrasts and fTests have the same number of columns as the
% design matrix
if ~isfield(params,'contrasts') || isempty(params.contrasts)
  params.contrasts = {};
else
  if size(d.scm,2)~=size(params.contrasts,2)*d.nHrfComponents
    mrErrorDlg( 'contrasts incompatible with number of EVs');
    return;
  end
end
contrasts = num2cell(params.contrasts,2); %convert contrasts to the same format as restrictions

if~isfield(params,'numberContrasts')
  params.numberContrasts = 0;
end
if~isfield(params,'computeTtests')
  params.computeTtests = 0;
end

if ~isfield(params,'tTestSide')
   params.tTestSide = 'Both';
end

if ~isfield(params,'fTestNames') || isempty(params.fTestNames)
   params.fTestNames = {};
end
  
if ~isfield(params,'restrictions') || isempty(params.restrictions)
   params.restrictions = {};
else
  if size(d.scm,2)~=size(params.restrictions{1},2)*d.nHrfComponents
    mrErrorDlg( 'F tests incompatible with number of EVs');
    return;
  end
  for iR = 1:length(params.restrictions)
    if ~all(any(params.restrictions{iR},2)) %remove empty contrasts from restriction matrix (not necessary but why not)
      params.restrictions{iR}(~any(params.restrictions{iR},2),:)=[];
    end
  end
end
restrictions = params.restrictions;

if ~isfield(d,'nHrfComponents')
  d.nHrfComponents = size(d.scm,2)/d.nhdr;
end
if ~isfield(params,'componentsToTest') || isempty(params.componentsToTest)
   params.componentsToTest = ones(1,d.nHrfComponents);
end

if ~isfield(params,'componentsCombination') || isempty(params.componentsCombination)
   params.componentsCombination = 'Add';
end
if strcmp(params.componentsCombination,'Or') && nnz(params.componentsToTest)==1
  params.componentsCombination = 'Add';  %make sure this is set to add if there is only one component to test
end

if fieldIsNotDefined(params,'outputStatistic')
  params.outputStatistic = 0;
end
params.outputStatistic = params.outputStatistic && actualData;

if fieldIsNotDefined(params,'bootstrapFweAdjustment')
  params.bootstrapFweAdjustment=0;
end
bootstrapFweAdjustment = params.bootstrapFweAdjustment && actualData;
if fieldIsNotDefined(params,'bootstrapTests')
  params.bootstrapTests=0;
end
if fieldIsNotDefined(params,'noDoubleBootstrap')
  params.noDoubleBootstrap=0;
end
bootstrapTests=params.bootstrapTests&& actualData;
if ~isfield(params, 'bootstrapIntervals') || isempty(params.bootstrapIntervals)
  params.bootstrapIntervals = 0;
end
bootstrapIntervals = params.bootstrapIntervals && actualData;
if ~isfield(params,'alphaConfidenceIntervals') || isempty(params.alphaConfidenceIntervals)
  params.alphaConfidenceIntervals = .05;
end

computeBootstrap = bootstrapTests || bootstrapFweAdjustment||bootstrapIntervals;
if ~computeBootstrap
  nResamples = 0;
elseif fieldIsNotDefined(params,'nResamples') 
  nResamples =10000;
else
  nResamples = params.nResamples;
end
if ieNotDefined('params') || ~isfield(params,'covCorrection') || ~isfield(params,'covCorrectionMethod') || strcmp(params.covCorrectionMethod,'none')
   params.covCorrection = 0;
end
if params.covCorrection
   if  ~isfield(params,'covCorrectionMethod') || isempty(params.covCorrectionMethod)
      params.covCorrectionMethod = 'preWhitening';
   end
   if  ~isfield(params,'covEstimation') || isempty(params.covEstimation)
      params.covEstimation = 'singleTukeyTapers';
   else
      params.covEstimation = params.covEstimation;
   end
   if  ~isfield(params,'covEstimationAreaSize') || isempty(params.covEstimationAreaSize)
      params.covEstimationAreaSize= 1;
   else
      params.covEstimationAreaSize = params.covEstimationAreaSize;
   end
   if  ~isfield(params,'covFactorization') || isempty(params.covFactorization)
      params.covFactorization= 'Cholesky';
   else
      params.covFactorization = params.covFactorization;
   end
else
   params.covCorrectionMethod = 'none';
end

%--------------------------------------------------- INITIALIZATION ------------------------------------%
%remove empty components if any
scm=d.scm;
scm(:,d.emptyEVcomponents)=[];
% precalculate the normal equation for Ordinary Least Squares
[invCovEVs,pinv_X]=computeNormalEquations(scm);

nFrames = size(d.volumes,2);

residualForming = eye(nFrames) - scm*pinv_X; %this matrix computes the residuals directly
d.rdf = nFrames-size(scm,2)-1; %degrees of freedom for the residual          %SHOULD BE nFrames-size(scm,2)

if ((params.covCorrection && params.covEstimationAreaSize>1) || params.spatialSmoothing)  && isfield(d,'roiPositionInBox') 
  d.dim(1) = nnz(d.roiPositionInBox);
end

%initialize variables for covariance matrix estimation
if params.covCorrection 
  if params.covEstimationAreaSize>1
    sliceAverager4D = ones(params.covEstimationAreaSize,params.covEstimationAreaSize)/params.covEstimationAreaSize^2;
    switch(params.covEstimationPlane)
      case {'Sagittal'}
        sliceAverager4D = permute(sliceAverager4D,[4 3 1 2]);
      case {'Axial'}
        sliceAverager4D = permute(sliceAverager4D,[4 1 2 3]);
      case {'Coronal'}
        sliceAverager4D = permute(sliceAverager4D,[4 1 3 2]);
    end
    if ~isfield(d,'roiPositionInBox') %this is in case we only compute for voxels in the loaded ROIs
      longMargin = ceil(params.covEstimationAreaSize/2); %we exclude voxels that are not surrounded by voxels with valid data
      shortMargin = floor(params.covEstimationAreaSize/2);
    end
  end
  switch params.covEstimation
    case 'singleTukeyTapers'
      tukeyM = 10;
      tukeyM = min(tukeyM,d.dim(4)-1);
      autoCorrelationParameters = NaN(tukeyM-1,d.dim(1),d.dim(2),d.dim(3),precision); %JB: this is to store the parameters of the ACF function (first 10 coeffs in this case)
  end
end

numberContrasts = length(contrasts);
if strcmp(params.componentsCombination,'Or') && nnz(params.componentsToTest)>1 %in this case, components of a given EV are tested independently from each other in an F-test
  params.componentsToTest = logical(diag(params.componentsToTest));
  %for contrasts, this amounts to  testing several contrats at once using an f test
  %So we have to change contrasts into f-tests (with the restriction that they have to be two-sided; this should be controlled for by the parameters)
  restrictions = [contrasts; restrictions];
  contrasts = {};
end
numberFtests = length(restrictions);

%expand restriction matrices using kronecker products
for iR = 1:numberFtests
  restrictions{iR} = kron(restrictions{iR},params.componentsToTest);
  %remove empty EV components
  restrictions{iR}(:,d.emptyEVcomponents)=[];
end
for iContrast = 1:length(contrasts)
  contrasts{iContrast} = kron(contrasts{iContrast},params.componentsToTest);
  %remove empty EV components
  contrasts{iContrast}(:,d.emptyEVcomponents)=[];
end


if ~isempty(restrictions)
  
  %this has not been tested, but this is only for generalized F-tests which do not give the expected results anyway
  if strcmp(params.covCorrectionMethod,'generalizedFTest')
    complementaryRestriction = cell(1,size(restrictions,1));
    baseRestriction =  kron(logical(eye(d.nhdr)),logical(params.componentsToTest)); 
    for iR = 1:numberFtests
      % not sure at all about these lines
      complementaryRestriction{iR} = baseRestriction - restrictions{iR};
      complementaryRestriction{iR} = complementaryRestriction{iR}(any(complementaryRestriction{iR},2),:);  
      % it actually doesn't make sense anymore now that F-tests can be made of any contrast
      % The logic relied on the fact that there were only elements on the diagonal of the restrictions matrix
      % There must be a generalization but I don't have time for this now
    end
    residualForming_f_tests = NaN(d.dim(4),d.dim(4),numberFtests);
    residualForming_h = NaN(d.dim(4),d.dim(4),numberFtests);
    traceRsV = NaN([d.dim(1:3) numberFtests]);
    effective_rdf = NaN(d.dim(1:3));
    effective_mdf = NaN([d.dim(1:3) numberFtests]);
    if ~approximate_value
      traceRV = NaN(d.dim(1:3));
    end
    scm_h = cell(1,numberFtests);
    eig_residualForming_f_tests = cell(1,numberFtests);
    for iR = 1:numberFtests
      scm_h{iR} = scm*complementaryRestriction{iR}';      
      residualForming_h(:,:,iR) = eye(nFrames) - scm_h{iR}*((scm_h{iR}'*scm_h{iR})^-1)*scm_h{iR}';       %TO REMOVE EVENTUALLY
      residualForming_f_tests(:,:,iR) = residualForming_h(:,:,iR) - residualForming; 
      V = eig(residualForming_f_tests(:,:,iR));
      eig_residualForming_f_tests{iR} = V(:,1:size(restrictions{iR},1));
    end
  end
  
  R_invCovEV_Rp = cell(1,numberFtests);
  for iR = 1:numberFtests         %the degrees of freedom for the F-tests are the number of contrasts
    restrictions{iR} = restrictions{iR}(any(restrictions{iR},2),:); %remove lines of 0
    d.mdf(iR) = size(restrictions{iR},1); %this is provided that contrasts are independent and that they are no more than (regressors -1)
    %inv_R_invCovEV_Rp{iR} = (restrictions{iR}*invCovEVs*restrictions{iR}')^-1; 
    %I replaced all the precomputed inverses (previous line) by ml/mrDivide with the non-inverted matrix (faster and more accurate):
    R_invCovEV_Rp{iR} = restrictions{iR}*invCovEVs*restrictions{iR}';
  end
                         
end


yvals = 1:d.dim(2);
slices = 1:d.dim(3);

numberTtests = length(contrasts);
numberTests = numberTtests+numberFtests;
  

%--------------------------------------------------- PRE-ALLOCATION ------------------------------------%

dInfo = whos('d');
totalBytes = dInfo.bytes; %keep track of big variables, including the raw data

[betas,totalBytes] = alloc('NaN',d.nhdr*d.nHrfComponents-length(d.emptyEVcomponents),d.dim(1),precision,totalBytes);
[thisS2,totalBytes] = alloc('NaN',d.dim(1),1,precision,totalBytes);
[thisResiduals,totalBytes] = alloc('NaN',d.dim(4),d.dim(1),precision,totalBytes);
if actualData 
  if bootstrapIntervals
    [sortedBetas,totalBytes] = alloc('NaN',d.nhdr*d.nHrfComponents-length(d.emptyEVcomponents),d.dim(1),nResamples,precision,totalBytes);    %JB: replaces: d.ehdr = zeros(d.dim(1),d.dim(2),d.dim(3),d.nhdr,d.hdrlen);
    [d.ehdrBootstrapCIs,totalBytes] = alloc('NaN',d.nhdr*d.nHrfComponents-length(d.emptyEVcomponents),d.dim(1),d.dim(2),d.dim(3),precision,totalBytes);
  end
  [d.ehdr,totalBytes] = alloc('NaN',d.nhdr*d.nHrfComponents-length(d.emptyEVcomponents),d.dim(1),d.dim(2),d.dim(3),precision,totalBytes);    %JB: replaces: d.ehdr = zeros(d.dim(1),d.dim(2),d.dim(3),d.nhdr,d.hdrlen);
  [rss,totalBytes] = alloc('NaN',d.dim(1),d.dim(2),d.dim(3),precision,totalBytes); %JB: this is to store the sum of square of the residual error term
  [tss,totalBytes] = alloc('NaN',d.dim(1),d.dim(2),d.dim(3),precision,totalBytes); %JB: this is to store the total sum of square
  [d.s2,totalBytes] = alloc('NaN',d.dim(1),d.dim(2),d.dim(3),precision,totalBytes); %JB: this is to store the estimated variance
end
if numberTtests
  [thisContrastBetas,totalBytes] = alloc('NaN',[d.dim(1) numberTtests],precision,totalBytes);
  [out.contrast,totalBytes] = alloc('NaN',[numberTtests d.dim(1) d.dim(2) d.dim(3)],precision,totalBytes);
  if bootstrapIntervals
    [sortedContrasts,totalBytes] = alloc('NaN',[numberTtests d.dim(1) nResamples],precision,totalBytes);
    [d.contrastBootstrapCIs,totalBytes] = alloc('NaN',[numberTtests d.dim(1) d.dim(2) d.dim(3)],precision,totalBytes);
  end
end
if numberFtests || params.computeTtests 
  [thisStatistic,totalBytes] = alloc('NaN',d.dim(1),numberTests,precision,totalBytes);
  [thisParametricP,totalBytes] = alloc('NaN',d.dim(1),numberTests,precision,totalBytes);
  [thisContrastBetaSte,totalBytes] = alloc('NaN',[d.dim(1) length(contrasts)],precision,totalBytes); 
  if computeBootstrap || params.outputStatistic
    [out.statistic,totalBytes] = alloc('NaN',[numberTests d.dim(1) d.dim(2) d.dim(3)],precision,totalBytes);
  end
  if (numberFtests ||params.computeTtests) && actualData
    %always output the parametric P if it's the actual data (even if params.parametricTests=0)
    %because it's fast and we might need it in some cases, but I don't want to do all the tests
    [out.parametricP,totalBytes] = alloc('NaN',[numberTests d.dim(1) d.dim(2) d.dim(3)],precision,totalBytes);
    if params.bootstrapFweAdjustment
      [out.bootstrapFweParametricP,totalBytes] = alloc('NaN',[numberTests d.dim(1) d.dim(2) d.dim(3)],precision,totalBytes);
    end
    if params.TFCE
      [out.bootstrapFweTfceP,totalBytes] = alloc('NaN',[numberTests d.dim(1) d.dim(2) d.dim(3)],precision,totalBytes);
    end
  end
  if bootstrapTests
    [out.bootstrapP,totalBytes] = alloc('NaN',[numberTests d.dim(1) d.dim(2) d.dim(3)],precision,totalBytes);
    if bootstrapFweAdjustment && ~params.noDoubleBootstrap
      [out.bootstrapFweBootstrapP,totalBytes] = alloc('NaN',[numberTests d.dim(1) d.dim(2) d.dim(3)],precision,totalBytes);
      [bootstrapStatistic,totalBytes] = alloc('NaN',[d.dim(1) numberTests nResamples],precision,totalBytes);
    end
    if params.TFCE
      [out.tfceBootstrapP,totalBytes] = alloc('NaN',[numberTests d.dim(1) d.dim(2) d.dim(3)],precision,totalBytes);
      if bootstrapFweAdjustment && ~params.noDoubleBootstrap
        [out.bootstrapFweTfceBootstrapP,totalBytes] = alloc('NaN',[numberTests d.dim(1) d.dim(2) d.dim(3)],precision,totalBytes);
        [bootstrapTfce,totalBytes] = alloc('NaN',[d.dim(1) numberTests nResamples],precision,totalBytes);
      end
    end
  end
  if numberFtests
    [mss,totalBytes] = alloc('NaN',d.dim(1),numberFtests,precision,totalBytes);
  end
  if params.covCorrection
    corrected_R_invCovEV_Rp = cell(1,numberFtests);
    for iR = 1:numberFtests
      [corrected_R_invCovEV_Rp{iR},totalBytes] = alloc('NaN',d.mdf(iR),d.mdf(iR),d.dim(1),'double',totalBytes);
    end
    switch(params.covCorrectionMethod)
      case 'generalizedLeastSquares'
        %pre allocate temporary variables
        [corrected_pinv_X,totalBytes] = alloc('NaN',d.nhdr*d.nHrfComponents-length(d.emptyEVcomponents),d.dim(4),d.dim(1),'double',totalBytes); %this array might become a bit large if a lot of values on the first dimension
        if actualData || numberTests
          [invCorrectedCovEV,totalBytes] = alloc('NaN',d.nhdr*d.nHrfComponents-length(d.emptyEVcomponents),d.nhdr*d.nHrfComponents-length(d.emptyEVcomponents),d.dim(1),'double',totalBytes);
        end
      case 'varianceCorrection' %see Woolrich et al. (2001) NeuroImage, 14(6), p1370
        [invCorrectedCovEV,totalBytes] = alloc('NaN',d.nhdr*d.nHrfComponents-length(d.emptyEVcomponents),d.nhdr*d.nHrfComponents-length(d.emptyEVcomponents),d.dim(1),'double',totalBytes);
        [correctedRdf,totalBytes] = alloc('NaN',d.dim(1),1,'double',totalBytes);
      case 'preWhitening' %see Woolrich et al. (2001) NeuroImage, 14(6), p1370
        [white_pinv_X,totalBytes] = alloc('NaN',d.nhdr*d.nHrfComponents-length(d.emptyEVcomponents),d.dim(4),d.dim(1),'double',totalBytes); %this array might become a bit large if a lot of voxels on the first dimension
        [white_scm,totalBytes] = alloc('NaN',d.dim(4),d.nhdr*d.nHrfComponents-length(d.emptyEVcomponents),d.dim(1),'double',totalBytes); %this array might become a bit large if a lot of voxels on the first dimension
        [correctedRdf,totalBytes] = alloc('NaN',d.dim(1),1,'double',totalBytes);
    end
  end
end

% turn off divide by zero warning
warning('off','MATLAB:divideByZero');



%--------------------------------------------------- MAIN LOOP: PARAMETER ESTIMATION ------------------------------------%
% display string
if verbose,hWaitBar = mrWaitBar(-inf,'(getGlmStatistics) Estimating model parameters');end
negativeSteWarning=false;
negativeS2Warning=false;
passesCounter = 0;
passesInBootstrap = d.dim(1);
switch(params.covCorrectionMethod)
  case 'none'
    preBootstrapTime = 1;
  case 'generalizedLeastSquares'
    preBootstrapTime = 1.4;
  case 'varianceCorrection'
    preBootstrapTime = d.dim(4)/2;
  case 'preWhitening'
    preBootstrapTime = d.dim(4);
end
totalPasses = prod(d.dim(1:3))*(nResamples+1+preBootstrapTime);

%these are the indices we will use to permute the residuals (with replacement)
if computeBootstrap %nResamples columns of randomly resampled indices 
  boostrapIndices = ceil(d.dim(4)*rand(nResamples,d.dim(4)));
end

% cycle through slices
for z = slices
  %permute timeseries frames in first dimension to compute residuals
  timeseries = permute(d.data(:,:,z,d.volumes),[4 1 2 3]); 
  % subtract off column means
  colmeans = repmat(mean(timeseries,1),[nFrames 1 1]);
  timeseries = timeseries - colmeans;
  % convert to percent signal change
  timeseries = 100*timeseries./colmeans;
  clear('colmeans');
  
  %spatial smoothing
  if params.spatialSmoothing
    if isfield(d,'roiPositionInBox') %if the data are not spatially organized, we need to temporarily put them in a volume
      timeseries = reshapeToRoiBox(timeseries,d.roiPositionInBox|d.marginVoxels,precision);
    end
    switch params.smoothingPlane %planes other than sagittal will only work for ROIs because it's the only case in which the data is 3D at this point in the loop
      case {'Sagittal'}
        timeseries = convn(timeseries,permute(gaussianKernel2D(params.spatialSmoothing),[4 3 1 2]),'same');
      case {'Axial'}
        timeseries = convn(timeseries,permute(gaussianKernel2D(params.spatialSmoothing),[4 1 2 3]),'same');
      case {'Coronal'}
        timeseries = convn(timeseries,permute(gaussianKernel2D(params.spatialSmoothing),[4 1 3 2]),'same');
      case '3D'
        timeseries = convn(timeseries,permute(gaussianKernel(params.spatialSmoothing),[4 1 2 3]),'same');
    end
    if isfield(d,'roiPositionInBox') %if the data are not spatially organized
      %put the data back in a new matrix with only the voxels of interest
      timeseries = reshape(timeseries,d.dim(4), numel(d.roiPositionInBox));
      if params.covCorrection
        timeseries = timeseries(:,d.roiPositionInBox|d.marginVoxels); %if we estimate the residuals covariance matrix, we keep the margin voxels
      else
        timeseries = timeseries(:,d.roiPositionInBox);
      end
      clear('usedVoxelsInBox');
    else
      
    end
  end    

  %--------------------------------------------------- MAIN LOOP: ORDINARY LEAST SQUARES ------------------------------------%
  % the following section has been optimized to run faster by
  % eliminating the x loops. in addition, the x dimension is divided in chunks of 60 points.
  % (which makes a difference only when system swaps ?)
  residuals = NaN(size(timeseries),precision);
%  residuals = NaN(d.dim(4),d.dim(1),d.dim(2),precision);
  for y = yvals
    % get OLS residuals
%     for iX = 1:ceil(d.dim(1)/optimalXLength)
%        xSubset = (iX-1)*optimalXLength+1 : min(iX*optimalXLength,d.dim(1));
%        residuals(:,xSubset,y) = residualForming*timeseries(:,xSubset,y); %JB resplaces residuals = timeseries-scm*d.ehdr(:,:,y,z);
%     end
    residuals(:,:,y) = residualForming*timeseries(:,:,y); %This is too long if X dim is large ?? (when calling with all data on the first dimension)
  end
  
  if params.covCorrection
    mrWaitBar( passesCounter/totalPasses, hWaitBar,'(getGlmStatistics) Estimating noise covariance (this can take a long time)');
    if params.covEstimationAreaSize>1     
      %average the residuals for covariance matrix estimation
      if isfield(d,'roiPositionInBox') %if the data are not spatially organized, we need to temporarily put them in a volume
%         %first put them on a single line whose length is the product of the volume dimensions
%         nVoxelsInBox = numel(d.roiPositionInBox);
%         usedVoxelsInBox = reshape(d.roiPositionInBox|d.marginVoxels,nVoxelsInBox,1);
%         averaged_residuals = NaN(d.dim(4),nVoxelsInBox,precision); 
%         averaged_residuals(:,usedVoxelsInBox) = residuals;  
%         %then reshape into a volume
%         averaged_residuals = reshape(averaged_residuals,[d.dim(4) size(d.roiPositionInBox)]);
        
        
        [averaged_residuals,usedVoxelsInBox] = reshapeToRoiBox(residuals,d.roiPositionInBox|d.marginVoxels,precision);
        %average
        if params.covCorrection && ~strcmp(params.covEstimationBrainMask,'None') && ~isempty(d.covEstimationBrainMask)
          averaged_residuals(:,~d.covEstimationBrainMask)=NaN;
          averaged_residuals = nanconvn(averaged_residuals,sliceAverager4D,'same'); %only do this if required because nanconvn takes a bit longer than convn
        else
          averaged_residuals = convn(averaged_residuals,sliceAverager4D,'same');
        end
        %put the data back in a new matrix with only the voxels of interest
        averaged_residuals = reshape(averaged_residuals,d.dim(4), numel(d.roiPositionInBox));
        averaged_residuals = averaged_residuals(:,d.roiPositionInBox);
        
        %reshape timeseries and residuals accordingly
        voxelsOfInterest = d.roiPositionInBox(usedVoxelsInBox);
        residuals = residuals(:,voxelsOfInterest);
        timeseries = timeseries(:,voxelsOfInterest);
        
        clear('usedVoxelsInBox');
      else
        averaged_residuals = NaN(size(residuals),precision);
        %averaged_residuals(isnan(residuals)) = 0;            %convert NaNs to 0... actually no, don't
        %average residuals over space
        averaged_residuals(:,longMargin:end-shortMargin,longMargin:end-shortMargin) = convn(residuals,sliceAverager4D,'valid');
      end
    else
      averaged_residuals = residuals;
    end
  end

  %---------------------- LOOP over Y: STATISTICS, NOISE COVARIANCE CORRECTION, RESIDUAL BOOTSTRAP ------------------------------------%
  for y = yvals  
 
    % Covariance matrix estimation 
    if params.covCorrection
      switch(params.covEstimation)
        %Single tapers
        case 'singleTukeyTapers'      %see Woolrich et al. (2001) NeuroImage, 14(6), p1370
          this_s2 = var(averaged_residuals(:,:,y),1,1);   
          for tau=1:tukeyM-1   
            autoCorrelationParameters(tau,:,y,z) = .5*(1+cos(pi*tau/tukeyM)) * (sum(averaged_residuals(1:end-tau,:,y).*averaged_residuals(tau+1:end,:,y),1))./this_s2/(d.dim(4)-tau);
          end
        case 'nonParametricResiduals'
            %not implemented
        case 'nonParametricTimeSeries' %Can we do this on the timeseries instead of the residuals, like Wicker et al ?
            %not implemented
        case 'dampenedOscillator' %see Kruggel et al. (2002) Medical image Analysis, 6, p65
           %check if invertible/definite positive (might have negative values ?)
           this_s2 = var(averaged_residuals(:,:,y),1,1);  
           M = d.dim(4)-1;
           %M=30;
           for tau=1:M
              autoCorrelation(tau+1,:) = (sum(averaged_residuals(1:end-tau,:,y).*averaged_residuals(tau+1:end,:,y),1))./this_s2/(d.dim(4)-tau);
           end
           %initialParameters = [0 0 0];
           initialParameters = rand(1,3);
           options = optimset('Display','iter');
           for x = xvals
              [dampenedOscillatorParams(x,:), dummy, exitflag] = fminsearch(@(dummy) minimizeDampenedOscillator(dummy,(1:M)',autoCorrelation(1+(1:M),x)), initialParameters,options);
           end
           %not finished ...
       end
    end
    
    if ~strcmp(params.covCorrectionMethod, 'generalizedFTest')
      %%%%%%%%%%%% COMPUTATIONS THAT CAN BE DONE BEFORE BOOTSTRAP (do not depend on the recomputed betas or residuals)
      if ismember(params.covCorrectionMethod,{'none','varianceCorrection'})
          %estimate the OLS beta weights 
          betas = pinv_X*timeseries(:,:,y); 
      end

      if ~strcmp(params.covCorrectionMethod, 'none')
        %find non-NaN voxels in this y line
        xvals = find(~any(isnan(timeseries(:,:,y))) & ~any(isnan(autoCorrelationParameters(:,:,y,z))));
        thisPassesCounter=passesCounter;
        for x = xvals
          if verbose,mrWaitBar( passesCounter/totalPasses, hWaitBar,'(getGlmStatistics) Estimating model parameters');end
          residualsAcm = makeAcm(autoCorrelationParameters(:,x,y,z),d.dim(4),params.covEstimation);
          
          switch(params.covCorrectionMethod)
          
            case 'generalizedLeastSquares' %see Wicker & Fonlupt (2003) NeuroImage, 18, p589 and Burock & Dale (2000) Human Brain Mapping, 11, p249
              correctedCovEV = scm' / residualsAcm * scm;
              %corrected_pinv_X(:,:,x) = correctedCovEV \ scm';
              %betas(:,x) = corrected_pinv_X(:,:,x)  / residualsAcm * timeseries(:,x,y);
              corrected_pinv_X(:,:,x) = correctedCovEV \ scm' / residualsAcm;
              betas(:,x) = corrected_pinv_X(:,:,x)   * timeseries(:,x,y);
              %these are the new GLS residuals (but we don't need the OLS ones anymore, so just replace)
              residuals(:,x,y) = timeseries(:,x,y) - scm*betas(:,x);
              if numberFtests || (numberTtests && params.computeTtests)
                invCorrectedCovEV(:,:,x) = inv(correctedCovEV); 
              end

            case 'varianceCorrection' %see Woolrich et al. (2001) NeuroImage, 14(6), p1370
              if numberFtests || (numberTtests&&params.computeTtests)
                invCorrectedCovEV(:,:,x) = pinv_X * residualsAcm * pinv_X'; 
              end
              correctedRdf(x) = trace(residualForming * residualsAcm);

            case 'preWhitening' %see Woolrich et al. (2001) NeuroImage, 14(6), p1370
              invResidualsAcm = inv(residualsAcm);
              % if this doesn't work, do pinv
              if sum(isnan(invResidualsAcm(:))) == length(invResidualsAcm(:))
                invResidualsAcm = pinv(residualsAcm);
                oneTimeWarning('pseudoInverseWarning','(getGlmStatistics) Residuals auto-correlation matrix is nearly singular, using pinv');
              end
              %compute inv(Cov)^1/2
              switch(params.covFactorization)
                case 'Cholesky'
                  try
                    preFilter = chol(invResidualsAcm);
                  catch exception
                    mrDisp(sprintf('(geGlmStatistics) Cannot factorize inverse noise covariance matrix because it is not positive definite (voxel %d %d %d)\n',x,y,z));
                    preFilter = [];
                  end
              end
              if ~isempty(preFilter)
                white_scm(:,:,x) = preFilter * scm;                 %%%% IS IT POSSIBLE TO SKIP THE FACTORIZATION BY ONLY USING (COV)^-1 LATER ON ??
                timeseries(:,x,y) = preFilter * timeseries(:,x,y);    %%%% NOT ACCORDING TO Kruggel et al. 2002...
                white_pinv_X(:,:,x) = (white_scm(:,:,x)'*white_scm(:,:,x))\white_scm(:,:,x)';
                white_residualForming = eye(nFrames) - white_scm(:,:,x)*white_pinv_X(:,:,x);
                betas(:,x) = white_pinv_X(:,:,x) * timeseries(:,x,y);
                if numberFtests || (numberTtests&&params.computeTtests) || nResamples>1
                  residuals(:,x,y) = white_residualForming * timeseries(:,x,y);
                end
                if numberFtests || (numberTtests&&params.computeTtests)
                  correctedRdf(x) = trace(white_residualForming * preFilter * residualsAcm * preFilter');
                  invCorrectedCovEV(:,:,x) = inv(scm' * invResidualsAcm * scm);
                end
              end
          end
            
          for iR = 1:numberFtests
            corrected_R_invCovEV_Rp{iR}(:,:,x) = restrictions{iR}*invCorrectedCovEV(:,:,x)*restrictions{iR}';
          end
          passesCounter = thisPassesCounter+x*preBootstrapTime;
        end
        passesCounter = thisPassesCounter+passesInBootstrap*preBootstrapTime;
      end
      
      %%%%%%%%%%%%%%%%% BOOTSTRAP LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for iBoot = 1:nResamples+1
        thisPassesCounter=passesCounter;
        if verbose 
          if iBoot==1
            mrWaitBar( passesCounter/totalPasses, hWaitBar,'(getGlmStatistics) Estimating model parameters');
          else
            mrWaitBar( passesCounter/totalPasses, hWaitBar,'(getGlmStatistics) Estimating model parameters (bootstrap)');
          end
        end
        if iBoot==1 
          thisResiduals = residuals(:,:,y);
          
          %initialize bootstrap counters
          if bootstrapTests
            bootstrapStatisticCount = zeros(d.dim(1), numberTests,precision);
            if params.TFCE
              bootstrapTfceCount = bootstrapStatisticCount;
            end
          end
          if (numberFtests ||params.computeTtests) && actualData && bootstrapFweAdjustment
            bootstrapFweStatisticCount = zeros(d.dim(1), numberTests,precision);
            if params.TFCE
              bootstrapFweTfceCount = bootstrapTfceCount;
            end
          end
        else
          % the bootstrap time series is resampled from the residuals with replacement. 
          % it SHOULD NOT be recomputed as beta*scm+residuals, as this violates pivotality 
          %(see Westfall, P.H., and S.S. Young. Resampling-based multiple testing. Wiley-Interscience, 1993. p35-39, p79-80, p108)
          % Basically, pivotality means that the statistics should be centered around 0 under the null hypothesis
          % Adding the GLM to the bootstrapped time-series would make the bootstrapped statistics 
          % centered around their estimate from the actual time-series. 
          % resulting in low power because they will be very close to the actual statistic value
          % On the other hand, if the GLM is included in actual time-series, not in the bootstrapped
          %  - if the null hypothesis is true for this test, both the actual statistic
          %        and the bootrapped statistics distribution will be close to 0
          %  - but if the null hypothesis is not true for this test, the actual statistic will be far
          %        from the bootstrapped statistic distribution (which will still be around 0)
          %
          % here we replace the timeseries by the bootstrap timeseries for this y, 
          % as we won't need the actual timeseries anymore
          timeseries(:,:,y) = residuals(boostrapIndices(iBoot-1,:),:,y);
          switch(params.covCorrectionMethod)
            case {'none','varianceCorrection'}
              %estimate the OLS beta weights 
              betas = pinv_X*timeseries(:,:,y); 
              %and we replace the residuals for this y by the bootstrap residuals
              thisResiduals = timeseries(:,:,y)-scm*betas;

            case 'generalizedLeastSquares' %see Wicker & Fonlupt (2003) NeuroImage, 18, p589 and Burock & Dale (2000) Human Brain Mapping, 11, p249
              for x = xvals
                if verbose,mrWaitBar( passesCounter/totalPasses, hWaitBar);end
                betas(:,x) = corrected_pinv_X(:,:,x) * timeseries(:,x,y);
                thisResiduals(:,x) = timeseries(:,x,y) - scm*betas(:,x);
                passesCounter = thisPassesCounter+x/2;
              end
              passesCounter = thisPassesCounter+passesInBootstrap/2;

            case 'preWhitening' %see Woolrich et al. (2001) NeuroImage, 14(6), p1370
              for x = xvals
                if verbose,mrWaitBar( passesCounter/totalPasses, hWaitBar);end
                betas(:,x) = white_pinv_X(:,:,x) * timeseries(:,x,y);
                thisResiduals(:,x) = timeseries(:,x,y) -white_scm(:,:,x)*betas(:,x);
                passesCounter = thisPassesCounter+x/2;
              end
              passesCounter = thisPassesCounter+passesInBootstrap/2;
          end
        end

        if actualData || numberFtests || (numberTtests&&params.computeTtests)
          switch(params.covCorrectionMethod)
            case 'none'
              thisRss = sum(thisResiduals.^2,1)';
              thisS2 = thisRss/d.rdf;
            case 'generalizedLeastSquares' %see Wicker & Fonlupt (2003) NeuroImage, 18, p589 and Burock & Dale (2000) Human Brain Mapping, 11, p249
              thisPassesCounter=passesCounter;
              for x = xvals
                if verbose,mrWaitBar( passesCounter/totalPasses, hWaitBar);end
                residualsAcm = makeAcm(autoCorrelationParameters(:,x,y,z),d.dim(4),params.covEstimation);
                %this is the part that makes GLS much slower than PW when bootstrapping
                thisS2(x) = thisResiduals(:,x)' / residualsAcm * thisResiduals(:,x)/d.rdf; 
%DEBUG                
%                 [~,messageId]=lastwarn;
%                 if strcmp(messageId,'MATLAB:nearlySingularMatrix')
%                   keyboard
%                 end
                 
                %Wicker & Fonlupt (2003)
                if thisS2(x)<0
                  thisS2(x)=NaN;
                  negativeS2Warning=true;
                end
                passesCounter = thisPassesCounter+x/2;
              end
              passesCounter = thisPassesCounter+passesInBootstrap/2;
            case {'varianceCorrection','preWhitening'} %see Woolrich et al. (2001) NeuroImage, 14(6), p1370
              thisRss = sum(thisResiduals.^2,1)';
              thisS2 = thisRss./correctedRdf;
          end
        end

        %----------------------------------  T-TESTS --------------------------------------------------------%
        if numberTtests 
          for iContrast = 1:numberTtests
            thisContrastBetas(:,iContrast) = (contrasts{iContrast}*betas)';
          end
          if params.computeTtests
            for iContrast = 1:numberTtests
              switch(params.covCorrectionMethod)
                case 'none'
                  thisContrastBetaSte(:,iContrast) = thisS2.*(contrasts{iContrast}*invCovEVs*contrasts{iContrast}');
                case {'generalizedLeastSquares','varianceCorrection','preWhitening'} 
                  for x = xvals
                    thisContrastBetaSte(x,iContrast) = thisS2(x)*contrasts{iContrast} * invCorrectedCovEV(:,:,x) * contrasts{iContrast}';
                    if thisContrastBetaSte(x,iContrast)<0
                      thisContrastBetaSte(x,iContrast)=NaN;
                      negativeSteWarning=true;
                    end
                  end
              end
            end
            thisContrastBetaSte = sqrt(thisContrastBetaSte);
            thisStatistic(:,1:numberTtests) = thisContrastBetas ./ thisContrastBetaSte;   
            switch(params.tTestSide)
             case 'Both'
                thisStatistic(:,1:numberTtests) = abs(thisStatistic(:,1:numberTtests));
             case 'Left'
                thisStatistic(:,1:numberTtests) = -1 *thisStatistic(:,1:numberTtests);
            end
            thisParametricP(:,1:numberTtests) = T2p(thisStatistic(:,1:numberTtests),d.rdf,params);
          end
        end
        %---------------------------------- F-TESTS --------------------------------------------------------%
        % basically, the F value (only for the EVs of interest) is 
        %     ((MSS)/ number of betas of interest) / (RSS/(number of samples - number of betas ?-1?))
        %     where MSS is the difference between the RSS without the EVs of interest and the RSS of the whole model
        % MSS has been computed as betas'*R'*(R*(X'*X)^-1*R')^-1 * R * betas, where R is a restrictions matrix that isolates the EVs of interest (such that R*X = Xinterest)
        % computed this way, the f-test can also be extended to contrasts and is equivalent to a two-sided T-test, if R describes a linear combination of EVs
        for iR = 1:numberFtests
          ss_beta = restrictions{iR}*betas;
          switch(params.covCorrectionMethod)
            case 'none'
              %computing time seems to be the fastest on my machine when computing 60 values at a time
              %but this might be because it was swapping. anyway, I'll leave it this way because if it doesn't swap
              %it doesn't seem to make much difference
      %                optimalXLength = 10:200;                             %DEBUG/OPTIMIZATION
      %                computingTimes = zeros(size(optimalXLength));        %DEBUG/OPTIMIZATION
      %                for j = 1:length(optimalXLength)                     %DEBUG/OPTIMIZATION
      %                   tic
              for iX = 1:ceil(d.dim(1)/optimalXLength)
                xSubset = (iX-1)*optimalXLength+1 : min(iX*optimalXLength,d.dim(1));
                %mss(xSubset) = diag( ss_beta(:,xSubset)' * inv_R_invCovEV_Rp{iR} * ss_beta(:,xSubset) ); 
                mss(xSubset,iR) = diag( ss_beta(:,xSubset)' / R_invCovEV_Rp{iR} * ss_beta(:,xSubset) );  
              end
      %                   computingTimes(j) = toc;                          %DEBUG/OPTIMIZATION
      %                end                                                  %DEBUG/OPTIMIZATION
      %                figure;plot(optimalXLength,computingTimes);          %DEBUG/OPTIMIZATION
              %mss = diag( ss_beta' * inv_R_invCovEV_Rp{iR} * ss_beta ); %%%%%this is too slow when Xdim is large
              %mss = diag( ss_beta' / R_invCovEV_Rp{iR} * ss_beta ); 

            case {'generalizedLeastSquares','varianceCorrection','preWhitening'} 
              for x = xvals
                 mss(x,iR) = ss_beta(:,x)' / corrected_R_invCovEV_Rp{iR}(:,:,x) * ss_beta(:,x); 
              end
          end
          thisStatistic(:,numberTtests+iR) = (mss(:,iR)/d.mdf(iR)) ./ thisS2;
        end
        if numberFtests
          thisParametricP(:,numberTtests+1:end) = F2p(thisStatistic(:,numberTtests+1:end),d.rdf,d.mdf);
        end
        
        if params.TFCE 
          thisTfce = applyTfce(thisStatistic,d.roiPositionInBox,precision); 
        end
        if numberFtests || params.computeTtests
          if iBoot==1
            actualStatistic = thisStatistic;
            if actualData && bootstrapFweAdjustment 
              resampleFweData = initResampleFWE(actualStatistic,params,thisParametricP);
            end
            if params.TFCE
              actualTfce = thisTfce; 
              if  actualData && bootstrapFweAdjustment
                resampleFweTfceTdata = initResampleFWE(actualTfce,params);
              end
            end
          else
            if bootstrapTests
              bootstrapStatisticCount = bootstrapStatisticCount+double(thisStatistic>actualStatistic); %T has already been transformed to be positive
              if bootstrapFweAdjustment && ~params.noDoubleBootstrap
                bootstrapStatistic(:,:,iBoot-1) = thisStatistic;
              end
              if params.TFCE
                bootstrapTfceCount = bootstrapTfceCount + double(thisTfce>actualTfce);
                if bootstrapFweAdjustment && ~params.noDoubleBootstrap
                  bootstrapTfce(:,:,iBoot-1) = thisTfce;
                end
              end
            end
            if actualData && bootstrapFweAdjustment
              bootstrapFweStatisticCount = bootstrapFweStatisticCount + resampleFWE(thisStatistic,actualStatistic,resampleFweData);
              if params.TFCE
                bootstrapFweTfceCount = bootstrapFweTfceCount + resampleFWE(thisTfce,actualTfce,resampleFweTfceTdata);
              end
            end
          end
        end

        if iBoot==1  %Actual values to save
          if actualData
            d.ehdr(:,:,y,z) = betas;
            switch(params.covCorrectionMethod)
              case 'none'
                rss(:,y,z) = thisRss;     
              case {'generalizedLeastSquares','varianceCorrection','preWhitening'} 
                rss(:,y,z) = sum(thisResiduals.^2,1)'; %SHOULD RSS BE CORRECTED ?
            end
            tss(:,y,z) = sum(timeseries(:,:,y).^2,1)';     
            d.s2(:,y,z) = thisS2;
          end
          if numberTtests
            out.contrast(:,:,y,z) = thisContrastBetas';
          end
          if (numberFtests ||params.computeTtests) && actualData
            out.parametricP(:,:,y,z) = thisParametricP';
          end
          if (numberFtests || (numberTtests&&params.computeTtests)) && params.outputStatistic
            out.statistic(:,:,y,z) = thisStatistic';
          end
        else
          if computeBootstrap && params.bootstrapIntervals && actualData
            % we keep all the beta estimates if we need to compute bootstrap confidence intervals
            sortedBetas(:,:,iBoot-1) = betas;
            if ~isempty(contrasts)
              sortedContrasts(:,:,iBoot-1) = thisContrastBetas';
            end
          end
        end
        
        if ismember(params.covCorrectionMethod,{'none','varianceCorrection','preWhitening'})
          passesCounter = thisPassesCounter+passesInBootstrap;
        end
      end
     
      %COMPUTE BOOTSTRAP STATISTICS
      if iBoot == nResamples+1 
        if numberTests
          if bootstrapTests
            bootstrapStatisticCount(isnan(actualStatistic))=NaN;
            out.bootstrapP(:,:,y,z) = computeBootstrapP(bootstrapStatisticCount',nResamples);
            if params.TFCE 
              bootstrapTfceCount(isnan(actualStatistic))=NaN;
              out.tfceBootstrapP(:,:,y,z) = computeBootstrapP(bootstrapTfceCount',nResamples);
            end
          end
          if actualData && bootstrapFweAdjustment
            bootstrapFweStatisticCount(isnan(actualStatistic))=NaN;
            bootstrapFweStatisticCount = computeBootstrapP(bootstrapFweStatisticCount,nResamples);
            out.bootstrapFweParametricP(:,:,y,z) = enforceMonotonicityResampleFWE(bootstrapFweStatisticCount,resampleFweData)';
            if params.TFCE 
              bootstrapFweTfceCount(isnan(actualStatistic))=NaN;
              bootstrapFweTfceCount = computeBootstrapP(bootstrapFweTfceCount,nResamples);
              out.bootstrapFweTfceP(:,:,y,z) = enforceMonotonicityResampleFWE(bootstrapFweTfceCount,resampleFweTfceTdata)';
            end
          end
          if bootstrapTests && bootstrapFweAdjustment && ~params.noDoubleBootstrap
           %SECOND BOOTSTRAP LOOP %one test at a time
            %convert to p-values
      % % %       [bootstrapStatistic, bootstrapTp] = sort(bootstrapStatistic,3);
      % % %       clear('bootstrapStatistic');
            [~, bootstrapStatistic] = sort(bootstrapStatistic,3); %this is not gonna work with older version of Matlab but is probably more efficient if bootstrapStatistic is really large
            %find order of actual values and estimate number of true null hypotheses
            resampleFweBootstrapData = initResampleFWE(out.bootstrapP(:,:,y,z)',params,out.bootstrapP(:,:,y,z)');
            out.bootstrapFweBootstrapP(:,:,y,z) = doubleResampleFWE(bootstrapStatistic,out.bootstrapP(:,:,y,z)',resampleFweBootstrapData,nResamples)';
            clear('bootstrapStatistic');
            out.bootstrapFweBootstrapP(:,:,y,z) = computeBootstrapP(out.bootstrapFweBootstrapP(:,:,y,z),nResamples);
            if params.TFCE
              %convert to p-values
              [~, bootstrapTfce] = sort(bootstrapTfce,3); %this is not gonna work with older version of Matlab but is probably more efficient if bootstrapStatistic is really large
              resampleFweBootstrapData = initResampleFWE(out.tfceBootstrapP(:,:,y,z)',params,out.tfceBootstrapP(:,:,y,z)');
              out.bootstrapFweTfceBootstrapP(:,:,y,z) = doubleResampleFWE(bootstrapTfce,out.tfceBootstrapP(:,:,y,z)',resampleFweBootstrapData,nResamples)';
              clear('bootstrapTfce');
              out.bootstrapFweTfceBootstrapP(:,:,y,z) = computeBootstrapP(out.bootstrapFweTfceBootstrapP(:,:,y,z),nResamples);
            end
          end
        end
        if params.bootstrapIntervals 
          lowerBoundIndex = max(1,floor(params.alphaConfidenceIntervals/2*nResamples));
          upperBoundIndex = min(nResamples,ceil((1-params.alphaConfidenceIntervals/2)*nResamples));
          sortedBetas = sort(sortedBetas,3);
          d.ehdrBootstrapCIs(:,:,y,z) = sortedBetas(:,:,upperBoundIndex)-sortedBetas(:,:,lowerBoundIndex);
          if numberTtests
            sortedContrasts = sort(sortedContrasts,3);
            d.contrastBootstrapCIs(:,:,y,z) = sortedContrasts(:,:,upperBoundIndex)-sortedContrasts(:,:,lowerBoundIndex);
          end
        end
      end
      
        
    else
      %Generalized F-tests - Not debugged
      %see Kruggel et al. (2002) Medical image Analysis, 6, p65
      %Does not give the expected result, but didn't have time to debug
      %anyway, less accurate and not obvious that it is much faster than other methods, so not worth the trouble
      d.ehdr(:,:,y,z) = pinv_X*timeseries(:,:,y);    % get OLS hdr
      if numberFtests 
        for x = xvals
          if verbose,mrWaitBar( (((y-min(yvals)) + d.dim(2)*(z-min(slices)))*d.dim(1)+x-min(xvals)), prod(d.dim(1:3)), hWaitBar);end
          if ~any(isnan(timeseries(:,x,y))) && ~any(isnan(autoCorrelation(:,x)))
                       residualsAcm = makeAcm(autoCorrelationParameters(:,x,y,z),d.dim(4),params.covEstimation); 
            if approximate_value 
               effective_rdf(x,y,z) = d.dim(4)*(d.dim(4)-d.nhdr*d.nHrfComponents+length(d.emptyEVcomponents))/trace(residualsAcm*residualsAcm);
            else %this is the exact value (longer and not very different for large N, which is usually the case for fMRI data)
               RV = residualForming*residualsAcm;           
               traceRV(x,y,z) = trace(RV);       
               effective_rdf(x,y,z) = traceRV(x,y,z)^2/trace(RV*RV);  
            end

            for iR = 1:numberFtests
               if approximate_value
                  temp = eig_residualForming_f_tests{iR}'*residualsAcm*eig_residualForming_f_tests{iR};   %SEE IF FASTER WHEN LOOP INSTEAD OF TRACE
                  effective_mdf(x,y,z,1,iR) = trace(temp)^2 / trace(temp.^2);                                      %SEE IF FASTER WHEN LOOP INSTEAD OF TRACE
                  traceRsV(x,y,z,1,iR) = trace(residualForming_f_tests(:,:,iR)*residualsAcm);
               else
                  RsV = residualForming_f_tests(:,:,iR)*residualsAcm;  %%this is the exact value (longer)
                  traceRsV(x,y,z,1,iR) = trace(RsV);
                  effective_mdf(x,y,z,1,iR) = traceRsV(x,y,z,1,iR)^2/trace(RsV*RsV);    
               end
            end
          end
        end
      end
      %from Burock and Dale 2000 
      for iR = 1:numberFtests
        ss_beta = restrictions{iR}*d.ehdr(:,:,y,z,1);
        for iX = 1:ceil(d.dim(1)/optimalXLength)
          xSubset = (iX-1)*optimalXLength+1 : min(iX*optimalXLength,d.dim(1));
          %mss(xSubset,y,z,1,iR) = diag( ss_beta(:,xSubset)' * inv_R_invCovEV_Rp{iR} * ss_beta(:,xSubset) ); 
          mss(xSubset,y,z,1,iR) = diag( ss_beta(:,xSubset)' / R_invCovEV_Rp{iR} * ss_beta(:,xSubset) ); 
        end
      %mss(:,y,z,1,iR) = diag( ss_beta' * inv_R_invCovEV_Rp{iR} * ss_beta );%%%%%%%%%%%%%%%%%%%%%%%%%%THIS IS LIKELY TO BE SLOW WHEN X IS LARGE
      %mss(:,y,z,1,iR) = diag( ss_beta' / R_invCovEV_Rp{iR} * ss_beta ); 
      end
      if approximate_value
        F = repmat(effective_rdf.^2,[1 1 1 1 numberFtests]) .* traceRsV .* mss ./ effective_mdf.^2 ./ repmat(rss,[1 1 1 1 numberFtests]) / (d.rdf+1);
      else
        F = repmat(effective_rdf.^2,[1 1 1 1 numberFtests]) .* traceRsV .* mss ./ effective_mdf.^2 ./ repmat(rss,[1 1 1 1 numberFtests]) ./ repmat(traceRV,[1 1 1 1 numberFtests]);
      end
    end
  end
  clear('averaged_residuals','residuals');
  clear('timeseries');
end
% disp([num2str(optimalXLength) ' ' num2str(t) ';']);   %DEBUG

clear('residualForming');
oneTimeWarning('pseudoInverseWarning',0);

%Now in the case we computed contrasts on several components using option 'Or',
%we have to convert the appropriate F values into T values
if numberContrasts>numberTtests  && params.outputStatistic
  %T values are the square roots of the numberContrasts first F values 
  out.statistic(1:numberContrasts,:,:,:) = sqrt(out.statistic(1:params.numberContrasts,:,:,:));
  d.mdf=d.mdf(params.numberContrasts+1:end);
end

%reshape to the right dimensions
if numberTtests
  out.contrast = permute(out.contrast,[2 3 4 1]);
end
if actualData
  % calculate variance accounted for by the estimated hdr
  out.r2 = 1 - rss./tss;     %JB: replaces: r2{y,z} = (1-sumOfSquaresResidual./sum(timeseries.^2));
  if ~isempty(d.emptyEVcomponents) %replace empty EV parameter estimates by nans
    nonEmptyEVcomponents = setdiff(1:d.nHrfComponents*d.nhdr,d.emptyEVcomponents);
    ehdr = nan(d.nHrfComponents*d.nhdr, d.dim(1), d.dim(2), d.dim(3));
    ehdr(nonEmptyEVcomponents,:,:,:) = d.ehdr;
    d.ehdr=ehdr;
  end
  %reshape nhdr and nhdrlen on two dimensions, put bootstrap dimension last, and swap nHrfComponents and nhdr order
  d.ehdr = permute(reshape(d.ehdr, d.nHrfComponents, d.nhdr, d.dim(1), d.dim(2), d.dim(3)),[3 4 5 2 1]);    
  if params.covCorrection
    d.autoCorrelationParameters = permute(autoCorrelationParameters,[2 3 4 1]);
  end
  if bootstrapIntervals
    if ~isempty(d.emptyEVcomponents)
      ehdrBootstrapCIs = nan(d.nHrfComponents*d.nhdr, d.dim(1), d.dim(2), d.dim(3));
      ehdrBootstrapCIs(nonEmptyEVcomponents,:,:,:) = d.ehdrBootstrapCIs;
      d.ehdrBootstrapCIs=ehdrBootstrapCIs;
    end
    d.ehdrBootstrapCIs = permute(reshape(d.ehdrBootstrapCIs, d.nHrfComponents, d.nhdr, d.dim(1), d.dim(2), d.dim(3)),[3 4 5 2 1]);
    if ~isempty(contrasts)
      d.contrastBootstrapCIs = permute(d.contrastBootstrapCIs,[2 3 4 1]);
    end
  end

end
if numberTests
  if (numberFtests ||params.computeTtests) && actualData
    out.parametricP = permute(out.parametricP,[2 3 4 1]);
    if bootstrapFweAdjustment
      out.bootstrapFweParametricP = permute(out.bootstrapFweParametricP,[2 3 4 1]);
      if params.TFCE 
        out.bootstrapFweTfceP = permute(out.bootstrapFweTfceP,[2 3 4 1]);
      end
    end
  end
  if (numberFtests || params.computeTtests) && params.outputStatistic
    out.statistic = permute(out.statistic,[2 3 4 1]);
  end
  if bootstrapTests
    out.bootstrapP = permute(out.bootstrapP,[2 3 4 1]);
    if bootstrapFweAdjustment && ~params.noDoubleBootstrap
      out.bootstrapFweBootstrapP = permute(out.bootstrapFweBootstrapP,[2 3 4 1]);
    end
    if params.TFCE 
      out.tfceBootstrapP = permute(out.tfceBootstrapP,[2 3 4 1]);
      if bootstrapFweAdjustment && ~params.noDoubleBootstrap
        out.bootstrapFweTfceBootstrapP = permute(out.bootstrapFweTfceBootstrapP,[2 3 4 1]);
      end
    end
  end
end
if verbose,mrCloseDlg(hWaitBar);end
if negativeS2Warning
  mrWarnDlg('(getGlmStatistics) There were some negative noise variance estimates');
end
if negativeSteWarning
  mrWarnDlg('(getGlmStatistics) There were some negative constrast variance estimates');
end


function [array, totalBytes] = alloc(defaultValues,varargin)

precision = varargin{end-1};
totalBytes = varargin{end};
dimensions = cell2mat(varargin(1:end-2));
arraySize = prod(dimensions);
switch(precision)
  case 'single'
    totalBytes = totalBytes+arraySize*4;
  case 'double'
    totalBytes = totalBytes+arraySize*8;
end

if totalBytes> mrGetPref('maxBlocksize')
  if ~askuser(sprintf('(getGlmStatistics) This analysis requires %.2f Gb of data to be loaded in memory at once,\nwhich is more than your MemoryBlock preference. This presumably happens because you are using a large ROI.\nConsider running the analysis on a (sub)volume or a smaller ROI.\nAre you sure you want to proceed ?', totalBytes/1024^3));
    error('(getGlmStatistics) Analysis aborted');
  end
end

switch(lower(defaultValues))
  case 'nan'
    array = NaN(dimensions,precision);
  case 'ones'
    array = ones(dimensions,precision);
  case 'zeros'
    array = zeros(dimensions,precision);
end

function [outputData,usedVoxelsInBox] = reshapeToRoiBox(inputData,dataPosition,precision)

%inputData must be size ([A B], where B equals the number of non-zero values in dataPosition
%outputData will be size ([A size(dataPosition)])
nVoxelsInBox = numel(dataPosition);
usedVoxelsInBox = reshape(dataPosition,nVoxelsInBox,1);
outputData = NaN(size(inputData,1),nVoxelsInBox,precision); 
outputData(:,usedVoxelsInBox) = inputData;  
%then reshape into a volume
outputData = reshape(outputData,[size(inputData,1) size(dataPosition)]);


function p = T2p(T,rdf,params)
p = 1 - cdf('t', double(T), rdf); %here use doubles to deal with small Ps
if strcmp(params.tTestSide,'Both')
  p = 2*p;
end
%we do not allow probabilities of 0 and replace them by minP
%(this can occur because cdf cannot return values less than 1e-16)
p = max(p,1e-16);
p(isnan(T)) = NaN; %NaNs must remain NaNs (they became 1e-16 when using max)


function p = F2p(F,rdf,mdf)
p = 1 - cdf('f', double(F), repmat(mdf,size(F,1),1), repmat(rdf,size(F)));  
%we do not allow probabilities of 0 and replace them by minP
%(this can occur because cdf cannot return values less than 1e-16)
p = max(p,1e-16);
p(isnan(F)) = NaN; %NaNs must remain NaNs (they became 1e-16 when using max)


function p = computeBootstrapP(count,nResamples)
%we do not allow probabilities of 0 and replace them by minP
%(this can occur if no bootstrap resampling return a value as high as the actual value)
p = max(count/nResamples,1/(nResamples+1));
p(isnan(count)) = NaN; %NaNs must remain NaNs (they became 1e-16 when using max)

function tfceS = applyTfce(S,roiPositionInBox,precision)
  %reshape to volume to apply TFCE and then reshape back to one dimension
  tfceS = applyFslTFCE(permute(reshapeToRoiBox(S',roiPositionInBox,precision),[2 3 4 1]),'',0);
  tfceS = permute(tfceS,[4 1 2 3]);
  tfceS = tfceS(:,roiPositionInBox)';
  %put NaNs back
  tfceS(isnan(S)) = NaN;

%this function computes the sum of squared errors between the dampened oscillator
%model (for xdata) and the sample autocorrelation function (ydata)
function sse = minimizeDampenedOscillator(params, xdata,ydata)
  FittedCurve = params(1)^2 - exp(params(2) * xdata) .* cos(params(3)*xdata);
  ErrorVector = FittedCurve - ydata;
  sse = sum(ErrorVector.^2);
  
 
  
function kernel = gaussianKernel(FWHM)

sigma_d = FWHM/2.35482;
w = ceil(FWHM); %deals with resolutions that are not integer
%make the gaussian kernel large enough for FWHM
kernelDims = 2*[w w w]+1;
kernelCenter = ceil(kernelDims/2);
[X,Y,Z] = meshgrid(1:kernelDims(1),1:kernelDims(2),1:kernelDims(3));
kernel = exp(-((X-kernelCenter(1)).^2+(Y-kernelCenter(2)).^2+(Z-kernelCenter(3)).^2)/(2*sigma_d^2)); %Gaussian function
kernel = kernel./sum(kernel(:));


function kernel = gaussianKernel2D(FWHM)

sigma_d = FWHM/2.35482;
w = ceil(FWHM); %deals with resolutions that are not integer
%make the gaussian kernel large enough for FWHM
kernelDims = 2*[w w]+1;
kernelCenter = ceil(kernelDims/2);
[X,Y] = meshgrid(1:kernelDims(1),1:kernelDims(2));
kernel = exp(-((X-kernelCenter(1)).^2+(Y-kernelCenter(2)).^2)/(2*sigma_d^2)); %Gaussian function
kernel = kernel./sum(kernel(:));
 

