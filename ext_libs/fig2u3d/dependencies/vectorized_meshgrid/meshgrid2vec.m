function [q] = meshgrid2vec(xgv, ygv, zgv)
%MESHGRID2VEC   meshgrid matrices to matrix of column vectors
%
%   [q] = MESHGRIDVEC(xgv, ygv) takes the matrix of abscissas XGV and
%   ordinates YGV of meshgrid points as returned by MESHGRID and arranges
%   them as vectors comprising the columns of matrix Q.
%
%   [q] = MESHGRIDVEC(xgv, ygv, zgv) does the same for the 3D case.
%
% input
%   xgv = matrix of meshgrid points' abscissas
%       = [#(points / y axis) x #(points / x axis) ] (2D case) or
%       = [#(points / y axis) x #(points / x axis) x #(points / z axis) ]
%         (3D case)
%   ygv = matrix of meshgrid points' ordinates
%       = similar dimensions with xgv
%   zgv = matrix of meshgrid points' coordinates
%       = similar dimensions with xgv
%
% output
%   q = [#dim x #(meshgrid points) ]
%
% See also DOMAIN2VEC, VEC2MESHGRID, DOMAIN2MESHGRID, MESHGRID.
%
% File:      meshgrid2vec.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.01.14 - 
% Language:  MATLAB R2011b
% Purpose:   convert meshgrid matrices to matrix of column vectors
% Copyright: Ioannis Filippidis, 2012-

% check input
if nargin == 2
    ndim = 2;
elseif nargin == 3
    ndim = 3;
else
    error('works only for ndim = 2 or ndim = 3/')
end

if ndim == 2
    q = meshgridvec2d(xgv, ygv);
elseif ndim == 3
    q = meshgridvec3d(xgv, ygv, zgv);
end

function [q] = meshgridvec2d(x, y)
q = [x(:), y(:) ].';

function [q] = meshgridvec3d(x, y, z)
q = [x(:), y(:), z(:) ].';
