function [vertices, edges] = line_pieces(vertices_all, piece_size, dashratio)
%
% See also u3d_pre_line, u3d_pre_contourgroup.
%
% File:      line_pieces.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.10 - 
% Language:  MATLAB R2012a
% Purpose:   generate vertices and edge arrays for pieces of a line
% Copyright: Ioannis Filippidis, 2012-

if nargin < 3
    dashratio = 0;
end

v = cut_line_to_pieces(vertices_all, piece_size);

n = size(v, 2);
vertices = cell(1, n);
edges = cell(1, n);
for i=1:n
    curv = v{1, i};
    
    %% connect to next line
    if i < n
        nv = size(curv, 2);
        next_line_v = v{1, i+1};
        vertex1_in_next_line = next_line_v(:, 1);
        
        curv(:, nv) = vertex1_in_next_line;
    end
    
    %% dashed style ?
    if 0 < dashratio
        nv = floor(piece_size /2);
        nv = max(1, nv);
        curv = curv(:, 1:nv); % keep half
    end
    
    %% create line
    npnt = size(curv, 2);
    curedges = [1:(npnt-1); 2:npnt]-1;
    
    %% output
    vertices{1, i} = curv;
    edges{1, i} = curedges;
end
