function [ThetaInDegrees, Column] = vectheta(varargin)
% Calculate the theta angle between two vectors or column vectors in a matrix

theta = @(u,v) atan2d(norm(cross(u,v)),dot(u,v));

if nargin == 1 % Input is a 3 x N matrix
    if isvector(varargin{1})
        error('Input must be a matrix');
    end
    A = varargin{1};
    Column = nchoosek(1:size(A,2), 2);
    ThetaInDegrees = zeros(size(Column,1), 1);
    for i=1:size(Column,1)
        u = A(:, Column(i,1));
        v = A(:, Column(i,2));
        ThetaInDegrees(i) = theta(u, v);
    end
elseif nargin == 2 % Input are two vectors
    if ~isvector(varargin{1}) || ~ isvector(varargin{2})
        error('Input must be vectors');
    end
    u = varargin{1};
    v = varargin{2};
    Column = [1 1];
    ThetaInDegrees = theta(u, v);
end
