function [X, Y, Z] = domain2meshgrid(domain, resolution)
%DOMAIN2MESHGRID(domain, resolution)   generate meshgrid on parallelepiped
%   [X, Y] = DOMAIN2MESHGRID(domain, resolution) creates the matrices
%   X, Y definining a meshgrid covering the 2D rectangular domain
%   domain = [xmin, xmax, ymin, ymax] with resolution = [nx, ny] points
%   per each coordinate dimension.
%
%usage
%-----
%   [X, Y, Z] = DOMAIN2MESHGRID(domain, resolution) results into a
%   meshgrid over a 3D parallelepiped domain.
%
%input
%-----
% (2D Case)
%   domain = extremal values of parallelepiped
%          = [xmin, xmax, ymin, ymax]
%   resolution = # points /dimension
%              = [nx, ny]
%
% (3D Case)
%   domain = [xmin, xmax, ymin, ymax, zmin, zmax]
%   resolution = [nx, ny, nz]
%
%output
%------
% (2D case)
%   X = [ny x nx] matrix of grid point abscissas
%   Y = [ny x nx] matrix of grid point ordinates
%
% (3D Case)
%   X = [ny x nx x nz] matrix of grid point abscissas
%   Y = [ny x nx x nz] matrix of grid point ordinates
%   Z = [ny x nz x nz] matrix of grid point coordinates
%
%about
%-----
%2012.01.14 (c) Ioannis Filippidis, jfilippidis@gmail.com
%
%See also DOMAIN2VEC, VEC2MESHGRID, MESHGRID2VEC, MESHGRID.

%% check input
if size(domain, 1) ~= 1
    error('size(domain, 1) ~= 1')
end

ndim = size(domain, 2) /2;

% if ~isint(ndim)
%     error('Non-integer domain dimension.')
% end

res_ndim = size(resolution, 2);
if res_ndim > ndim
    warning('dom2meshgrid:res_ndim',...
          ['size(resolution) = [', num2str(size(resolution) ),...
           '] ~= [', num2str([1, ndim] ), '] = [1, ndim]'] )
end

if res_ndim < ndim
    error(['size(resolution) = [', num2str(size(resolution) ),...
           '] ~= [', num2str([1, ndim] ), '] = [1, ndim]'] )
end

if ndim == 2
    [X, Y] = linmeshgrid2d(domain, resolution);
elseif ndim == 3
    [X, Y, Z] = linmeshgrid3d(domain, resolution);
else
    msg = 'domain has more than 3 dimensions. Use vec2 directly.';
    warning('dom:dim4', msg)
end

function [X, Y] = linmeshgrid2d(domain, resolution)
xmin = domain(1, 1);
xmax = domain(1, 2);

ymin = domain(1, 3);
ymax = domain(1, 4);

nx = resolution(1, 1);
ny = resolution(1, 2);

x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);

[X, Y] = meshgrid(x, y);

function [X, Y, Z] = linmeshgrid3d(domain, resolution)
xmin = domain(1, 1);
xmax = domain(1, 2);

ymin = domain(1, 3);
ymax = domain(1, 4);

zmin = domain(1, 5);
zmax = domain(1, 6);

nx = resolution(1, 1);
ny = resolution(1, 2);
nz = resolution(1, 3);

x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);
z = linspace(zmin, zmax, nz);

[X, Y, Z] = meshgrid(x, y, z);
