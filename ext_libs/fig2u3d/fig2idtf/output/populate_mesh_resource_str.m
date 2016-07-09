function [mesh_resources] = populate_mesh_resource_str(f, v, c)
%
% See also FACE_VERTEX_DATA_EQUALS_NPOINTS, MESH_NORMALS, VERBATIM.
%
% File:      populate_mesh_resource_str.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.10 - 
% Language:  MATLAB R2012a
% Purpose:   mesh resource strings
% Copyright: Ioannis Filippidis, 2012-
%
% acknowledgment
%   Based on save_idtf by Alexandre Gramfort.
%   This can be found on the MATLAB Central File Exchange:
%       http://www.mathworks.com/matlabcentral/fileexchange/25383-matlab-mesh-to-pdf-with-3d-interactive-object
%   and is covered by the BSD License.

% depends
%   single_mesh_resource_str

n_meshes = size(f, 2);
tmp_mesh_resources = cell(1, n_meshes);
for i=1:n_meshes
    faces = f{1, i};
    points = v{1, i};
    face_vertex_data = c{1, i};
    
    cur_mesh_resource = single_mesh_resource_str(faces, points, face_vertex_data, i);
            
    tmp_mesh_resources{1, i} = ['\n\n', cur_mesh_resource];
end
mesh_resources = [tmp_mesh_resources{:} ];
%disp(mesh_resources)
