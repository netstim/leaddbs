% makeAcm.m
%
%        $Id: computeNormalEquations.m 1950 2010-12-18 10:12:48Z julien $	
%      usage: residualsAcm = makeAcm(autoCorrelationParameters,sampleSize,covEstimationType)
%         by: julien besle
%       date: 04/01/2011
%    purpose: computes residual autocorrelation matrix for timeseries of length sampleSize
%                       from parameters estimated using method covEstimationType  
%                       (only 'singleTukeyTapers' implemented for now)
%             
function residualsAcm = makeAcm(autoCorrelationParameters,sampleSize,covEstimationType)

switch(covEstimationType)
  case 'singleTukeyTapers'
    autoCorrelation = zeros(sampleSize,1,'double'); %you don't want to make those single
    %because inversion operations are faster on doubles
    autoCorrelation(1) = 1;
    autoCorrelation(2:size(autoCorrelationParameters,1)+1) = autoCorrelationParameters;
    residualsAcm = toeplitz(autoCorrelation);
end