function [d, V] = center(X, V, dim)
%CENTER subtracts a vector from each row of a matrix
%
%data utility to simplify subtracting each row of a matrix by a vector
%
% d = center(X, V, dim)
% X is a m x n matrix
% V is a 1 x n vector. defaults to column sum of X
% dim is a scalar equal to 1 or 2. 1 to subtract a row vector, 2
% to subtract a column vector
% d is a new matrix with each row of X divided by V
%
%Example
%   load carbig
%   X = [MPG Acceleration Weight Displacement];
%   z = center( X, nanmean(X) ); 
%   z = center(X, mean(X,2), 2);

% $Id: center.m,v 1.4 2006/12/26 22:54:06 Mike Exp $
% Copyright 2006 Mike Boedigheimer
% Amgen Inc.
% Department of Computational Biology
% mboedigh@amgen.com
% 

if nargin < 3
    dim = 1;
end;

if ( nargin < 2 )
    V = mean(X, dim);
end;


[m, n] = size(X);
s = n;
if dim == 2
    s = m;
end;

if ~isscalar(V) && (~isvector(V) || length(V) ~= s)
     error('linstats:scale:InvalidArgument', 'Input argument must be a vector');
end


if dim==2
    d = X - repmat(V, 1, n);
else
    d = X - repmat(V, m, 1);
end;