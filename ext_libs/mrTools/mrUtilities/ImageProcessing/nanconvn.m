% nanconvn - performs n-dimensional convolution, ignoring NaNs in A and weighting each point 
%            by the actual number of datapoints involved (weighted by B values)
%
%        $Id$
%      usage: [C,weights] = nanconvn(A,B,shape)
%         by: julien besle
%       date: 2011-02-10
%     inputs: same as convn
%       goal: values whose computation involve NaNs are computed, replacing NaNs with zeros
%             and normalizing according to the weighted number of non-NaN values involved
%
%             Note that for simplicity, values outside the valid area are also 
%             modified relative to what convn would output. This is because
%             less points are involved in the computation of the edges (because of zero padding in convn)
%             This could be corrected by padding the border of the mask with ones, 
%             but would require a bit more involvment in bookkeeping the valid/full shape
%             As not normalizing the edges (like in convn) seems as arbitrary 
%             as normalizing, I couldn't be bothered.
%
%             This function only accounts for NaNs in array A

function [C,weights] = nanconvn(A,B,shape)


if ieNotDefined('shape')
  shape = 'full';
end

%create a mask of non-NaN values
mask = ~isnan(A);
%convolving the mask with B will give weights representing how many and how much
%actual data points actually contributed to the convolved values 

%replace NaNs by zeros in A and B
A(isnan(A))=0;

%convolve thsi new A with B and divide by the weights
weights=convn(mask,abs(B),shape);
C = convn(A,B,shape)./weights;

%Note that some points will be undefined (those that don't receive any contribution from any point of A)



