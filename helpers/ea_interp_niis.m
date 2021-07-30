function niiinterp=ea_interp_niis(nii,interstages)
niiinterp=[];
for stage=1:length(nii)-1
    niiinterp=cat(4,niiinterp,linspaceNDim(nii(stage).img,nii(stage+1).img,interstages));
end



function y = linspaceNDim(d1, d2, n)
%LINSPACENDIM Linearly spaced multidimensional matrix.
%   LINSPACENDIM(d1, d2) generates a multi-dimensional 
%   matrix of 100 linearly equally spaced points between 
%   each element of matrices d1 and d2.
%
%   LINSPACENDIM(d1, d2, N) generates N points between 
%   each element of matrices X1 and X2.
%
%       Example:
%       d1 = rand(3, 2, 4); d2 = rand(size(d1)); n = 10; 
%
%       y = linspaceNDim(d1, d2, n) returns a multidimensional matrix y of
%       size (3, 2, 4, 10)
%
%
%   Class support for inputs X1,X2:
%      float: Multidimensional matrix, vector, double, single
%
%       Steeve AMBROISE --> sambroise@gmail.com
%
%   $ Date: 2009/01/29   21:00:00 GMT $ 
%   $ revised Date: 2009/02/02   18:00:00 GMT $ 
%       Bug fixed for singleton dimensions that occur when d1 or d2 
%           are empty matrix, scalar or vector.
% 
%
if nargin == 2
    n = 100;
end
n  = double(n);
d1 = squeeze(d1); d2 = squeeze(d2);
if ndims(d1)~= ndims(d2) || any(size(d1)~= size(d2))
    error('d1 and d2 must have the same number of dimension and the same size'),
end
NDim = ndims(d1);
%%%%%%%% To know if the two first dimensions are singleton dimensions
if NDim==2 && any(size(d1)==1)
    NDim = NDim-1;
    if all(size(d1)==1)
        NDim = 0;
    end
end
pp      = (0:n-2)./(floor(n)-1);
Sum1    = TensorProduct(d1, ones(1,n-1));
Sum2    = TensorProduct((d2-d1), pp);
y = cat(NDim+1, Sum1  + Sum2, shiftdim(d2, size(d1, 1)==1 ));
%%%%% An old function that I wrote to replace the built in Matlab function:
%%%%% KRON
function Z = TensorProduct(X,Y)
%   Z = TensorProduct(X,Y) returns the REAL Kronecker tensor product of X and Y. 
%   The result is a multidimensional array formed by taking all possible products
%   between the elements of X and those of Y. 
%
%   If X is m-by-n and Y is p-by-q-by-r, then kron(X,Y) 
%   is m-by-p-by-n-by-q-by-r.
%
%   X and Y are multidimensional matrices 
%   of any size and number of dimensions 
% 
%   E.g. if X is of dimensions (4, 5, 3) and Y of dimension (3, 1, 7, 4)
%   TensorProduct(X, Y) returns a multidimensional matrix Z of dimensions: 
%   (4, 5, 3, 3, 7, 4)
% 
%   $ Date: 2001/11/09 10:20:00 GMT $ 
%
%       Steeve AMBROISE --> sambroise@gmail.com
%
sX=size(X);sY=size(Y);
ndim1=ndims(X);ndim2=ndims(Y);
indperm=[ndim2+1:ndim1+ndim2,1:ndim2];
% to remove all singleton dimensions 
Z=squeeze(repmat(X,[ones(1,ndims(X)),sY]).*... 
    permute(repmat(Y,[ones(1,ndims(Y)),sX]),indperm));