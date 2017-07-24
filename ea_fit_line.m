function [trajectory,trajvector,estpoint]=ea_fit_line(centerline)
% This function fits a line to raw 3D data specified by centerline. It also
% exports a normalized vector between the fitted datapoints and the next
% expected point of the trajectory. The function uses robustmean by
% Developer "Jonas"
% (http://www.mathworks.com/matlabcentral/fileexchange/[...]
% 27134-plot-average-line/content/plotAverage/robustMean.m)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

    trajectory=centerline;
    points=[0:size(trajectory,1)-1]';

    linelen=size(centerline,1)-1;
    gaussweights=zeros(linelen,1);
    for i=1:linelen
        gaussweights(i)=exp(-((i-linelen/2)^2)/(0.7*linelen));
    end
    
    
    xversatz=ea_robustmean(diff(centerline(1:end,1))); %wmean(diff(centerline(1:end,1)),gaussweights,1);
    yversatz=ea_robustmean(diff(centerline(1:end,2)));%wmean(diff(centerline(1:end,2)),gaussweights,1);
    
    
    
    
    correctstartx=ea_robustmean(centerline(1:end,1)-[points].*xversatz);
    correctstarty=ea_robustmean(centerline(1:end,2)-[points].*yversatz);
    
    trajectory=[correctstartx+points.*xversatz,correctstarty+points.*yversatz,centerline(1,3)-points];
    

    
    trajvector=[xversatz,yversatz,-1];
    
    
    % estimated next point.
    estpoint=trajectory(end,:)+trajvector;
    
function y = ea_nanmean(varargin)
if nargin==2
    x=varargin{1};
    dim=varargin{2};
elseif nargin==1
x=varargin{1};
    dim=1;
end
    
N = sum(~isnan(x), dim);
y = sum(x(~isnan(x)), dim) ./ N;

   
    
   

function [finalMean, stdSample, inlierIdx, outlierIdx] = ea_robustmean(data,dim,k,fit)
%ROBUSTMEAN calculates mean and standard deviation discarding outliers
%
% SYNOPSIS [finalMean, stdSample, inlierIdx, outlierIdx] = robustMean(data,dim,k,fit)
%
% INPUT    data : input data
%          dim  : (opt) dimension along which the mean is taken {1}
%          k    : (opt) #of sigmas at which to place cut-off {3}
%          fit  : (opt) whether or not to use fitting to robustly estimate
%                  the mean from the data that includes outliers.
%                  0 (default): mean is approximated by median(data)
%                  1 : mean is approximated by
%                      fminsearch(@(x)(median(abs(data-x))),median(data))
%                      This option is only available for scalar data
%
%
% OUTPUT   finalMean : robust mean
%          stdSample : std of the data (divide by sqrt(n) to get std of the
%                      mean)
%          inlierIdx : index into data with the inliers
%          outlierIdx: index into data with the outliers
%
% REMARKS  NaN or Inf will be counted as neither in- nor outlier
%          The code is based on (linear)LeastMedianSquares. It could be changed to
%          include weights
%
% c: jonas, 04/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% test input

if isempty(data)
    error('Please supply non-empty data to robustMean')
end
if nargin<2 || isempty(dim)
    % make sure that the dimensinon is correct if there's a vector
    if any(size(data)==1) && ndims(data)==2
        dim = find(size(data)>1);
    else
        dim = 1;
    end
end
if nargin < 3 || isempty(k)
    k = 3;
end
if nargin < 5 || isempty(fit)
    fit = 0;
end
if fit == 1
    % check for vector
    if sum(size(data)>1)>1
        error('fitting is currently only supported for 1D data')
    end
end
if sum(isfinite(data(:))) < 4
    finalMean = ea_nanmean(data,dim);
    stdSample = NaN(size(finalMean));
    inlierIdx = find(isfinite(data));
    outlierIdx = [];
    return
end


%========================
% LEAST MEDIAN SQUARES
%========================

% define magic numbers:
%k=3; %cut-off is roughly at 3 sigma, see Danuser, 1992 or Rousseeuw & Leroy, 1987
magicNumber2=1.4826^2; %see same publications

% remember data size and reduced dataSize
dataSize = size(data);
reducedDataSize = dataSize;
reducedDataSize(dim) = 1;
% need this for later repmats
blowUpDataSize = dataSize./reducedDataSize;
% count how many relevant dimensions we have besides dim
realDimensions = length(find(dataSize>1));

% calc median - reduce dimension dim to length 1
if fit
    % minimize the median deviation from the mean
    medianData = fminsearch(@(x)(median(abs(data-x))),median(data));
else
    medianData = ea_nanmedian(data,dim);
end

% calculate statistics
res2 = (data-repmat(medianData,blowUpDataSize)).^2;
medRes2 = max(ea_nanmedian(res2,dim),eps);

%testvalue to calculate weights
testValue=res2./repmat(magicNumber2*medRes2,blowUpDataSize);

if realDimensions == 1;
    %goodRows: weight 1, badRows: weight 0
    inlierIdx=find(testValue<=k^2);
    outlierIdx = find(testValue>k^2);
    
    % calculate std of the sample;
    if nargout > 1
        nInlier = length(inlierIdx);
        if nInlier > 4
            stdSample=sqrt(sum(res2(inlierIdx))/(nInlier-4));
        else
            stdSample = NaN;
        end
    end
    
    %====END LMS=========
    
    %======
    % MEAN
    %======
    
    finalMean = mean(data(inlierIdx));
    
else
    
    %goodRows: weight 1, badRows: weight 0
    inlierIdx=find(testValue<=k^2);
    outlierIdx=find(testValue > k^2);
    
    % mask outliers
    res2(outlierIdx) = NaN;
    % count inliers
    nInliers = sum(~isnan(res2),dim);
    
    % calculate std of the sample;
    if nargout > 1
        % put NaN wherever there are not enough data points to calculate a
        % standard deviation
        goodIdx = sum(isfinite(res2),dim) > 4;
        stdSample = NaN(size(goodIdx));
        stdSample(goodIdx)=sqrt(nansum(res2(goodIdx),dim)./(nInliers(goodIdx)-4));
    end
    
    %====END LMS=========
    
    %======
    % MEAN
    %======
    
    data(outlierIdx) = NaN;
    finalMean = ea_nanmean(data,dim);
end
