function [q, X, Y, Z] = domain2vec(domain, resolution)
%DOMAIN2VEC     matrix of column vectors for meshgrid points
%   [q, X, Y, Z] = domain2vec(domain, resolution) returns a matrix of
%   column vectors for the meshgrid point coordinates over a parallelepiped
%   domain with the given resolution. The individual coordinates are also
%   returned in separate matrices X, Y, Z.
%
% input (2D Case)
%   domain = extremal values of parallelepiped
%          = [xmin, xmax, ymin, ymax]
%   resolution = # points /dimension
%              = [nx, ny]
%
% input (3D Case)
%   domain = [xmin, xmax, ymin, ymax, zmin, zmax]
%   resolution = [nx, ny, nz]
%
% output
%   q = matrix of column vectors (meshgrid point coordinates)
%     = [#dim x #points]
%   X = meshgrid point abscissas (nz = 1 for the 2D case)
%     = [ny x nx x nz]
%   Y = meshgrid point ordinates (nz =1 for the 2D case)
%     = [ny x nx x nz]
%   Z = meshgrid point coordinates (defined only for the 3D case)
%     = [ny x nx x nz]
%
% See also VEC2MESHGRID, DOMAIN2MESHGRID, MESHGRID2VEC.
%
% File:      domain2vec.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2010.09.16 - 2012.01.22
% Language:  MATLAB R2010b
% Purpose:   
% Copyright: Ioannis Filippidis, 2010-

ndim = size(domain, 2) /2;

if ndim == 2
    [X, Y] = domain2meshgrid(domain, resolution);
    q = meshgrid2vec(X, Y);
    if nargout == 4
        Z = zeros(size(X) );
        warning('domain2vec:Z', 'Z matrix is 0 for 2D grid points.')
    end
elseif ndim == 3
    [X, Y, Z] = domain2meshgrid(domain, resolution);
    q = meshgrid2vec(X, Y, Z);
else
    error('ndim \notin {2, 3}')
end
