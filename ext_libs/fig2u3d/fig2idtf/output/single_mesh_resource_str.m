function [mesh_resource] = single_mesh_resource_str(faces, points, face_vertex_data, i)
%
% See also u3d_pre_patch.
%
% File:      single_mesh_resource_str.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.10 - 2012.06.24
% Language:  MATLAB R2012a
% Purpose:   export single mesh to idtf
% Copyright: Ioannis Filippidis, 2012-
%
% acknowledgment
%   Based on save_idtf by Alexandre Gramfort.
%   This can be found on the MATLAB Central File Exchange:
%       http://www.mathworks.com/matlabcentral/fileexchange/25383-matlab-mesh-to-pdf-with-3d-interactive-object
%   and is covered by the BSD License.

% depends
%   mesh_diffuse_colors, mesh_normals, verbatim

points = points.';

nfaces = size(faces, 1);
npoints = size(points, 1);
nface_vertex_data = size(face_vertex_data, 1);
resource_number = i -1;

% RGB cdata ?
if ~isempty(face_vertex_data) && size(face_vertex_data, 2) ~= 3
    error('Colors should be RGB.');
end

% #cdata == #points ?
switch nface_vertex_data
    case 0
        disp('Empty cdata.')
    case npoints
        disp('size(cdata, 1) = #points')
    case nfaces
        disp('size(cdata, 1) = #faces')
    otherwise
        disp(['#colors = ', num2str(nface_vertex_data) ] )
        disp(['#points = ', num2str(npoints) ] )
        disp(['#faces = ', num2str(nfaces) ] )
        keyboard
        error('#colors ~= #points and #colors ~= #faces and #colors not empty.')
end

%% colors
[strmodelshading, strfaces_colors, strfaceshading, strcolors, nface_vertex_data_unique] = ...
    mesh_diffuse_colors(faces, points, face_vertex_data);

strfaces = sprintf('                    %d %d %d\n', faces.'-1);
strpoints = sprintf('                    %f %f %f\n', points.');

[normals, ~] = compute_normal(points, faces);
normals(isnan(normals) ) = 0;
strnormals = sprintf('                    %f %f %f\n', normals.');

% model diffuse color and shading counts
if (nface_vertex_data == npoints) || (nface_vertex_data == 0)
    model_diffuse_color_count = nface_vertex_data_unique;
    model_shading_count = 1;
elseif nface_vertex_data == nfaces
    model_diffuse_color_count = 0;
    model_shading_count = nface_vertex_data_unique;
else
    error('#colors ~= #points, #faces and not empty. (previously checked!)')
end

mesh_resource = sprintf(resource_mesh_str, resource_number, i,...
            nfaces, npoints, npoints,...
            model_diffuse_color_count, model_shading_count,...
            strmodelshading, strfaces, strfaces, strfaceshading,...
            strfaces_colors, strpoints, strnormals, strcolors);

function [str] = resource_mesh_str
% resource_number, mesh_number, face_count, model_position_count,
% model_normal_count,
% model_diffuse_color_count, model_shading_count,
% model_shading_description_list_contents
%
% mesh_face_position_list, mesh_face_normal_list,
% mesh_face_shading_list, mesh_face_diffuse_color_list,
%
% model_position_list, model_normal_list, model_diffuse_color_list
str = verbatim;
%{
    RESOURCE %d {
          RESOURCE_NAME "MyMesh%d"
          MODEL_TYPE "MESH"
          MESH {
               FACE_COUNT %d
               MODEL_POSITION_COUNT %d
               MODEL_NORMAL_COUNT %d
               MODEL_DIFFUSE_COLOR_COUNT %d
               MODEL_SPECULAR_COLOR_COUNT 0
               MODEL_TEXTURE_COORD_COUNT 0
               MODEL_BONE_COUNT 0
               MODEL_SHADING_COUNT %d
               MODEL_SHADING_DESCRIPTION_LIST {
%s
               }
               MESH_FACE_POSITION_LIST {
%s
               }
               MESH_FACE_NORMAL_LIST {
%s
               }
%s
               
%s
               MODEL_POSITION_LIST {
%s
               }
               MODEL_NORMAL_LIST {
%s
               }
%s
          }
     }
%}
