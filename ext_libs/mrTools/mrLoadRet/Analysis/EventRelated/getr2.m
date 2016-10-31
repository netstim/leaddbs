% getr2.m
%
%      usage: getr2(d)
%         by: justin gardner
%       date: 07/30/03
%    purpose: function to get r2 image by
%             finding hdr for each voxel
%             and seeing how much of the
%             variance the estimate accounts for.
%       e.g.: d = getr2(d)
%             d needs to have a stimulus convolution matrix
%             in d.scm
%
function d = getr2(d,verbose)

% init some variables
ehdr=[];r2 = [];

if ieNotDefined('verbose'),verbose = 1;end

% precalculate the normal equation (this dramatically speeds up things)
covarianceMatrix = d.scm'*d.scm;
covarianceMatrixRank = rank(covarianceMatrix);
% use normal equations if we have a full ranked covariance matrix
if covarianceMatrixRank == size(covarianceMatrix,1)
  inverseCovarianceMatrix = covarianceMatrix^-1;
  precalcmatrix = inverseCovarianceMatrix*d.scm';
  diagOfInverseCovariance = diag(inverseCovarianceMatrix);
% otherwise use pinv
else
  % note that if we need to use the pseudo inverse it means that there is ambiguity in the design
  % such that there are an infinite number of possible solutions. The psuedo-inverse solution
  % chosses the solution with the minimum length (i.e. Euclidian norm)
  if verbose,disp(sprintf('(getr2) Design covariance matrix (%ix%i) is rank %i. Using pseudo-inverse to invert.',size(covarianceMatrix,1),size(covarianceMatrix,2),covarianceMatrixRank));end
  precalcmatrix = pinv(d.scm);
  % get the diagonal of the inverse of the design covariance matrix (used for estimating standard errors)
  diagOfInverseCovariance = diag(pinv(d.scm'*d.scm));
end

% check roi
slices = 1:d.dim(3);slicen = length(slices);
xvals = 1:d.dim(1);xvaln = length(xvals);
yvals = 1:d.dim(2);yvaln = length(yvals);
  
% preallocate memory
d.ehdr = zeros(d.dim(1),d.dim(2),d.dim(3),d.nhdr,d.hdrlen);
d.ehdrste = zeros(d.dim(1),d.dim(2),d.dim(3),d.nhdr,d.hdrlen);
d.r2 = zeros(d.dim(1),d.dim(2),d.dim(3));

% check for image dimensions that this code will choke on
% this should not normally happen, but can happen if you
% are using maxBlocksize to limit the size of images that
% are brought in at a time
  if isequal(d.dim(1),1) && isequal(d.dim(3),1)
    mrWarnDlg(sprintf('(getr2) Dimensions of image are degenerate: %i x %i x %i (should have more than one dimension with size greater than 1. This may be caused by your maxBlocksize setting in Preferences. Try setting it larger to allow larger parts of your data to be loaded into memory at the same time',d.dim(1),d.dim(2),d.dim(3)));
    keyboard
  end

% turn off warnings to avoid divide by zero warning
warning('off','MATLAB:divideByZero');

% display string
if verbose,disppercent(-inf,'(getr2) Calculating r2');end
% cycle through images calculating the estimated hdr and r^2s of the 
% estimate.
%
% this following section has been optimized to run faster by
% eliminating one of the loops. Various different methods were
% tested eliminating all the loops and doing one big calculation
% which thrashed memory too much, or eliminating different
% dimensions and it was found that eliminating the first dimension
% was by far the faster by a factor of about 2-3. 
onesmatrix = ones(length(d.volumes),1);
for j = yvals
  for k = slices
    % get the time series we are working on
    % this includes all the rows of one column from one slice
    % and all data points for each of these
    % thus the time series is a nxm matrix where each of the m columns
    % contains the n time points recording for that voxel
    timeseries = squeeze(d.data(:,j,k,d.volumes))';
    % subtract off column means
    colmeans = mean(timeseries,1);
    timeseries = timeseries - onesmatrix*colmeans;
    % convert to percent signal change
    timeseries = 100*timeseries./(onesmatrix*colmeans);
    % get hdr for the each voxel
    ehdr{j,k} = precalcmatrix*timeseries;
    % calculate error bars, first get sum-of-squares of residual
    % (in percent signal change)
    sumOfSquaresResidual = sum((timeseries-d.scm*ehdr{j,k}).^2);
    % now calculate the sum-of-squares of that error
    % and divide by the degrees of freedom (n-k where n
    % is the number of timepoints in the scan and k is 
    % the number of timepoints in all the estimated hdr)
    S2 = sumOfSquaresResidual/(length(d.volumes)-size(d.scm,2));
    % now distribute that error to each one of the points
    % in the hemodynamic response according to the inverse
    % of the covariance of the stimulus convolution matrix.
    ehdrste{j,k} = sqrt(diagOfInverseCovariance*S2);
    % calculate variance accounted for by the estimated hdr
    r2{j,k} = (1-sumOfSquaresResidual./sum(timeseries.^2));
  end
  if verbose,disppercent(max((j-min(yvals))/yvaln,0.1));end
end
if verbose,disppercent(inf);end

% reshape matrix. this also seems the fastest way to do things. we
% could have made a matrix in the above code and then reshaped here
% but the reallocs needed to continually add space to the matrix
% seems to be slower than the loops needed here to reconstruct
% the matrix from the {} arrays.
if verbose,disppercent(-inf,'(getr2) Reshaping matrices');end
for i = xvals
  for j = yvals
    for k = slices
      % get the ehdr
      d.ehdr(i,j,k,1:d.nhdr,:) = reshape(squeeze(ehdr{j,k}(:,i)),[d.hdrlen d.nhdr])';
      % and the stderror of that
      d.ehdrste(i,j,k,1:d.nhdr,:) = reshape(squeeze(ehdrste{j,k}(:,i)),[d.hdrlen d.nhdr])';
      % now reshape r2 into a matrix
      d.r2(i,j,k) = r2{j,k}(i);
    end
  end
  if verbose,disppercent((i-min(xvals))/xvaln);end
end

% display time took
if verbose,disppercent(inf);end

warning('on','MATLAB:divideByZero');
