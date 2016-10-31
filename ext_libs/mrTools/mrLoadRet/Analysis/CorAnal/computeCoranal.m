% modelCoranalFromGlm
%
%        $Id$
%      usage: [co, amp, ph] = computeCoranal(tSeries,params)
%         by: julien besle
%       date: 2010-11-15
%    purpose: computes correlation analysis on timeseries array, called by corAnal, corAnalPlot and modelCoranalFromGLM
%        e.g:
%
function [co, amp, ph, ptSeries] = computeCoranal(tSeries,nCycles,detrend,spatialNormalization,trigonometricFunction)

nFrames = size(tSeries,1);
% Set highpassPeriod
highpassPeriod = round(nFrames/nCycles);

% Remove dc, convert to percent, detrend, and spatial normalization
warnState = warning('query','MATLAB:divideByZero');
ptSeries = percentTSeries(tSeries,...
    'detrend', detrend,...
    'highpassPeriod', highpassPeriod,...
    'spatialNormalization', spatialNormalization,...
    'subtractMean', 'Yes',...
    'temporalNormalization', 'No');
warning(warnState.state,warnState.identifier);

% Compute Fourier transform
ft = fft(ptSeries);
ft = ft(1:1+fix(size(ft, 1)/2), :);
ampFT = 2*abs(ft)/nFrames;

% Compute co and amp (avoiding divide by zero)
amp = ampFT(nCycles+1,:);
co = zeros(size(amp),'single');
sumAmp = sqrt(sum(ampFT.^2));
nonzeroIndices = find(sumAmp >0);
co(nonzeroIndices) = ampFT(nCycles+1,nonzeroIndices) ./ sumAmp(nonzeroIndices);

% Calculate phase:
switch(trigonometricFunction)
  case 'Sine'
    % 1) add pi/2 so that it is in sine phase.
    % 2) minus sign because sin(x-phi) is shifted to the right by phi.
    ph = -(pi/2) - angle(ft(nCycles+1,:));   %
  case 'Cosine'
    ph = - angle(ft(nCycles+1,:));   
end
% 3) Add 2pi to any negative values so phases increase from 0 to 2pi.
ph(ph<0) = ph(ph<0)+pi*2;