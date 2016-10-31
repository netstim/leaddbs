% getEstimates.m
%
%        $Id: getEstimates.m 1950 2010-12-18 10:12:48Z julien $	
%      usage: getEstimates(d,params,indices,contrasts)
%         by: julien besle
%       date: 03/01/2011
%    purpose: retrieve/compute the estimate, contrast estimates and their standard error from d structure. 
%                indices must be a n*3 matrix of [x y s] subscripts or a vector of linear indices indexing the volume 
%                (there is an ambiguity if indices is a 1*3 vector, in which case, it is considered a single vector of subscripts)
function [estimates,outputIndices] = getEstimates(d,params,indices)

if ~any(nargin == [2 3])
  help getContrastEstimate;
  return
end

% if x is a list of 3-tuples
if size(indices,2) == 3
  indices = sub2ind(d.dim(1:3),indices(:,1),indices(:,2),indices(:,3));
elseif length(size(indices))==2 && any(size(indices) == 1)
  if size(indices,1) == 1
    indices = indices';
  end
else
  mrWarnDlg('(getEstimates) Wrong voxel indices format');
  return
end
nVoxels = length(indices);

if fieldIsNotDefined(d,'estimationSupersampling')
  d.estimationSupersampling =1;
end
if fieldIsNotDefined(d,'designSupersampling')
  d.designSupersampling=1;
end
if fieldIsNotDefined(d,'acquisitionDelay')
  d.acquisitionDelay=d.tr/2;
end

betas = permute(d.ehdr,[4 5 1 2 3]);
betas = betas(:,:,indices);
noValue = isnan(squeeze(betas(1,1,:)));

if fieldIsNotDefined(d,'actualhrf')
  if fieldIsNotDefined(d,'hrf')
    hrf = eye(size(betas,2));
    params.componentsToTest = ones(1,size(betas,2));
  else
    hrf = d.hrf;
  end
else
  hrf = d.actualhrf;
end
hrfLength = size(hrf,1);

if fieldIsNotDefined(d,'nHrfComponents')
  d.nHrfComponents = size(hrf,2);
end

if isfield(params,'componentsCombination')
  nComponents = d.nHrfComponents*d.nhdr;
  nHrfComponents = d.nHrfComponents;
else
  nComponents = d.nhdr;
  nHrfComponents = 1;
end

extendedHrf = kron(eye(d.nhdr),hrf);
betas = reshape(permute(betas,[2 1 3]),nComponents,nVoxels);
if fieldIsNotDefined(d,'emptyEVcomponents')
  hdr = extendedHrf*betas;
else
  nonEmptyEVcomponents = setdiff(1:d.nHrfComponents*d.nhdr,d.emptyEVcomponents);
  extendedHrf(:,d.emptyEVcomponents)=[];
  hdr = extendedHrf*betas(nonEmptyEVcomponents,:);
end
betas = permute(reshape(betas,[nHrfComponents d.nhdr nVoxels]),[2 1 3]);
hdr = reshape(hdr,[hrfLength,d.nhdr,nVoxels]);

% Compute estimate std error (see for example Seber & Lee, Linear Regression Analysis, 1977, p42,67)
hdrSte = NaN(size(hdr));
if ~fieldIsNotDefined(d,'s2')
  %remove empty components if any
  if ~fieldIsNotDefined(d,'emptyEVcomponents')
    d.scm(:,d.emptyEVcomponents)=[];
  end
  [invCovEVs,pinv_X]=computeNormalEquations(d.scm);
  s2 = permute(d.s2(indices),[2 3 1]);
  %for hdrSte:
  if params.covCorrection && ~fieldIsNotDefined(d,'autoCorrelationParameters')
    betaSte = NaN(size(betas));
    acfParams = permute(d.autoCorrelationParameters,[4 1 2 3]);
    acfParams = acfParams(:,indices);
    invCorrectedCovEV = NaN(size(d.scm,2),size(d.scm,2),nVoxels);
    for iVoxel = 1:nVoxels
      if ~isnan(s2(iVoxel))
        residualsAcm = makeAcm(acfParams(:,iVoxel),d.dim(4),params.covEstimation);
        %THIS IS WHAT TAKES TIME
        switch(params.covCorrectionMethod)
          case {'generalizedLeastSquares','preWhitening'} 
            %GLS: see Wicker & Fonlupt (2003) NeuroImage, 18, p589 and Burock & Dale (2000) Human Brain Mapping, 11, p249
            %PW: see Woolrich et al. (2001) NeuroImage, 14(6), p1370
            invCorrectedCovEV(:,:,iVoxel) = inv(d.scm' / residualsAcm * d.scm); 
          case 'varianceCorrection' %see Woolrich et al. (2001) NeuroImage, 14(6), p1370
            invCorrectedCovEV(:,:,iVoxel) = pinv_X * residualsAcm * pinv_X'; 
        end
        if ~fieldIsNotDefined(d,'emptyEVcomponents')
          thisBetaSte = nan(nHrfComponents*d.nhdr,1);
          thisBetaSte(nonEmptyEVcomponents)=diag(invCorrectedCovEV(:,:,iVoxel));
        else
          thisBetaSte = diag(invCorrectedCovEV(:,:,iVoxel));
        end
        betaSte(:,:,iVoxel) = reshape(thisBetaSte,nHrfComponents,d.nhdr)';
        % according to Friston et al. 1998 Neuroimage 7 p30-40:
        hdrSte(:,:,iVoxel) = reshape(diag(extendedHrf*invCorrectedCovEV(:,:,iVoxel)*extendedHrf'),hrfLength,d.nhdr);
      end
    end
  else %OLS 
    if ~fieldIsNotDefined(d,'emptyEVcomponents')
      betaSte = nan(nHrfComponents*d.nhdr,1);
      betaSte(nonEmptyEVcomponents) = diag(invCovEVs);
    else
      betaSte = diag(invCovEVs);
    end
    betaSte = repmat(reshape(betaSte,nHrfComponents,d.nhdr)',[1 1 nVoxels]);
    hdrSte = repmat(reshape(diag(extendedHrf*invCovEVs*extendedHrf'),hrfLength,d.nhdr),[1 1 nVoxels]);
  end
  betaSte = sqrt(betaSte.*repmat(s2,[d.nhdr,nHrfComponents 1]));  
  hdrSte = sqrt(hdrSte.*repmat(s2,[hrfLength d.nhdr 1]));  
  s2 = permute(s2,[3 1 2]);
  noValue = isnan(s2);
elseif ~fieldIsNotDefined(d,'ehdrste') && size(betas,1)==size(d.ehdrste,4) && size(betas,2)==size(d.ehdrste,5)
  betaSte = permute(d.ehdrste,[4 5 1 2 3]);
  betaSte = betaSte(:,:,indices);
  for iVoxel = 1:nVoxels
    %assuming that components of the HRF are independent (which always is true for deconvolution, but might not for other models)
    %the variance (squared stddev) of a random variable multiplied by a constant 
    %is equal to the variance multiplied by the square of the constant. 
    hdrSte(:,:,iVoxel) = sqrt((hrf.^2)*(betaSte(:,:,iVoxel).^2)'); 
  end
end


estimates.contrastBetas = [];
estimates.contrastBetaSte = [];
estimates.contrastHdr = [];
estimates.contrastHdrSte = [];
if ~fieldIsNotDefined(params,'contrasts')
  if size(params.contrasts,2)~= size(betas,1)
     mrWarnDlg(sprintf('(getEstimate) The number of columns in the contrast vector should equal the number of EVs (%d)',size(betas,1)));
  elseif isfield(params,'componentsCombination') && length(params.componentsToTest)~=d.nHrfComponents
     mrWarnDlg(sprintf('(getEstimate) The number of EV components to test is misspecified'));
  else
    nContrasts = size(params.contrasts,1);
    if isfield(params,'componentsCombination')
      %estimated std error for hdr contrast  will be estimated separately
      hdrContrasts = kron(params.contrasts,hrf*diag(params.componentsToTest));
      %for the beta estimates, we use either the 'Add' or the 'Or' mode
      switch(params.componentsCombination)
        case 'Add'
          extendedContrasts = kron(params.contrasts,params.componentsToTest);
        case 'Or'
          extendedContrasts = kron(params.contrasts,diag(params.componentsToTest)); 
          %remove empty contrasts
          extendedContrasts = extendedContrasts(any(extendedContrasts,2),:);
      end
    else
      hdrContrasts = params.contrasts;
      extendedContrasts = params.contrasts;
    end

    betas = reshape(permute(betas,[2 1 3]),nComponents,nVoxels);
    if ~fieldIsNotDefined(d,'emptyEVcomponents')
      extendedContrasts = extendedContrasts(:,nonEmptyEVcomponents);
      hdrContrasts = hdrContrasts(:,nonEmptyEVcomponents);
      contrastBetas = extendedContrasts*betas(nonEmptyEVcomponents,:);
      contrastHdr = reshape(hdrContrasts*betas(nonEmptyEVcomponents,:),hrfLength,nContrasts,nVoxels);
    else
      contrastBetas = extendedContrasts*betas;
      contrastHdr = reshape(hdrContrasts*betas,hrfLength,nContrasts,nVoxels);
    end
    betas = permute(reshape(betas,[nHrfComponents d.nhdr nVoxels]),[2 1 3]);

    contrastBetaSte = NaN(size(contrastBetas));
    contrastHdrSte = NaN(size(contrastHdr));
    if exist('s2','var')
      if params.covCorrection  && ~fieldIsNotDefined(d,'autoCorrelationParameters')
        contrastHdrSte = NaN(size(contrastHdr));
        for iVoxel = 1:nVoxels
          contrastBetaSte(:,iVoxel) = sqrt(diag(extendedContrasts * invCorrectedCovEV(:,:,iVoxel) * extendedContrasts')*s2(iVoxel));
          contrastHdrSte(:,:,iVoxel) = sqrt(reshape(diag(hdrContrasts * invCorrectedCovEV(:,:,iVoxel) * hdrContrasts'),hrfLength,nContrasts)*s2(iVoxel));
        end
      else
        contrastBetaSte = sqrt(diag(extendedContrasts*invCovEVs*extendedContrasts')*s2');
        contrastHdrSte = sqrt(reshape(diag(hdrContrasts*invCovEVs*hdrContrasts')*s2',hrfLength,nContrasts,nVoxels));
  %      contrastHdrSte = reshape(contrastHdrSte,hrfLength,nContrasts,nVoxels);
      end

      contrastBetaSte = permute(reshape(contrastBetaSte,[],nContrasts,nVoxels),[2 1 3]);

    elseif ~fieldIsNotDefined(d,'contrastSte') && size(extendedContrasts,1)==size(d.contrastSte,4)*size(d.contrastSte,5)
      contrastBetaSte = permute(d.contrastSte,[4 5 1 2 3]);
      contrastBetaSte = contrastBetaSte(:,:,indices);
    end

    %reshape in case there are several contrast components
    contrastBetas = permute(reshape(contrastBetas,[],nContrasts,nVoxels),[2 1 3]);
    contrastBetaSte = permute(reshape(contrastBetaSte,[],nContrasts,nVoxels),[2 1 3]);
    
    %copy in output structure
    estimates.contrastBetas = contrastBetas(:,:,~noValue);
    estimates.contrastBetaSte = contrastBetaSte(:,:,~noValue);
    estimates.contrastHdr = contrastHdr(:,:,~noValue);
    estimates.contrastHdrSte = contrastHdrSte(:,:,~noValue);
  end

end

outputIndices = indices(~noValue);
estimates.betas = betas(:,:,~noValue);
estimates.betaSte = betaSte(:,:,~noValue);
%set any null estimate to NaN
hdr(hdr==0)=NaN;
hdrSte(hdrSte==0)=NaN;
estimates.hdr = hdr(:,:,~noValue);
estimates.hdrSte = hdrSte(:,:,~noValue);


%mean confidence intervals across voxels 
if isfield(d,'ehdrBootstrapCIs')
  estimates.bootstrapBetaCIs = permute(d.ehdrBootstrapCIs,[4 5 1 2 3]);
  estimates.bootstrapBetaCIs = estimates.bootstrapBetaCIs(:,:,indices(~noValue));
  if ~fieldIsNotDefined(d,'contrastBootstrapCIs') && isequal(d.contrasts,params.contrasts)
    estimates.bootstrapContrastCIs = permute(d.contrastBootstrapCIs,[4 5 1 2 3]);
    estimates.bootstrapContrastCIs = estimates.bootstrapContrastCIs(:,:,indices(~noValue));
  else
    estimates.bootstrapContrastCIs = NaN(size(estimates.contrastBetas));
  end
end

% Time vector
tr = d.tr/max(d.estimationSupersampling,d.designSupersampling);
% estimates.time = ((d.acquisitionSubsample-.5)*tr:tr:(hrfLength+d.acquisitionSubsample-1.5)*tr)';
estimates.time = rem(d.acquisitionDelay,tr)+tr*(0:hrfLength-1);

