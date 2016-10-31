% computeNormalEquations.m
%
%        $Id: computeNormalEquations.m 1950 2010-12-18 10:12:48Z julien $	
%      usage: [invCovEVs,pinv_X]=computeNormalEquations(designMatrix)
%         by: julien besle
%       date: 04/01/2011
%    purpose: 
%             

function [invCovEVs,pinv_X]=computeNormalEquations(designMatrix)

covEVs = designMatrix'*designMatrix; %covariance matrix of the EVs (do I need to correct for covariance ?)
covarianceMatrixRank = rank(covEVs);
% use normal equations if we have a full ranked covariance matrix
if covarianceMatrixRank == size(covEVs,1)
   invCovEVs = covEVs^-1; % 
   pinv_X = covEVs\designMatrix'; 
        
% otherwise use pinv
else
  % note that if we need to use the pseudo inverse it means that there is ambiguity in the design
  % such that there are an infinite number of possible solutions. The pseudo-inverse solution
  % choses the solution with the minimum length (i.e. Euclidian norm)
  mrWarnDlg(sprintf('(computeNormalEquations) Design covariance matrix (%ix%i) is rank %i. Using pseudo-inverse to invert.',size(covEVs,1),size(covEVs,2),covarianceMatrixRank))
  pinv_X = pinv(designMatrix);
  invCovEVs = pinv(designMatrix'*designMatrix);
end
