function [strmodelshading, strfaces_colors, strfaceshading, strcolors,...
          nface_vertex_data_unique] = mesh_diffuse_colors(faces, points, face_vertex_data)
% File:      mesh_diffuse_colors.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.10- 2012.06.24
% Language:  MATLAB R2012a
% Purpose:   prepare colors and shading for exporting mesh to IDTF
% Copyright: Ioannis Filippidis, 2012-
%
% acknowledgment
%   Based on save_idtf by Alexandre Gramfort.
%   This can be found on the MATLAB Central File Exchange:
%       http://www.mathworks.com/matlabcentral/fileexchange/25383-matlab-mesh-to-pdf-with-3d-interactive-object
%   and is covered by the BSD License.

% depends
%   verbatim

%% no colors
if isempty(face_vertex_data)
    strmodelshading = sprintf(shading_description_str, [0, 0].');
    strfaces_colors = '';
    strfaceshading = sprintf(mesh_face_shading_list, '0');
    strcolors = '';
    nface_vertex_data_unique = 0;
    return
end

%% colors
nfaces = size(faces, 1);
npoints = size(points, 1);
nfvd = size(face_vertex_data, 1);

[face_vertex_data_unique, ~, face_vertex_data_idx] = ...
                                unique(face_vertex_data, 'rows');
nface_vertex_data_unique = size(face_vertex_data_unique, 1);
msg = ['Number of colors: ', num2str(nface_vertex_data_unique) ];
disp(msg);

% Colors correspond to points or faces ?
if nfvd == npoints
    strmodelshading = sprintf(shading_description_str, [0, 0].');
    
    strfaces_colors = sprintf('                %d %d %d\n', face_vertex_data_idx(faces).'-1);
    strfaces_colors = sprintf(mesh_face_diffuse_color_list, strfaces_colors);
    
    strfaceshading = sprintf('                    %d\n', zeros(nfaces, 1) );
    strfaceshading = sprintf(mesh_face_shading_list, strfaceshading);
    
    strcolors = sprintf('                %f %f %f\n', face_vertex_data_unique.');
    strcolors = sprintf(model_diffuse_color_list, strcolors);
elseif nfvd == nfaces
    shidx = (0:(nface_vertex_data_unique -1) ).';
    strmodelshading = sprintf(shading_description_str, [shidx, shidx].');
    
    strfaces_colors = '';
    
    strfaceshading = sprintf('                %d\n', face_vertex_data_idx -1);
    strfaceshading = sprintf(mesh_face_shading_list, strfaceshading);
    
    strcolors = '';
else
    error('#facevetrexdata ~= #faces, #points, 0')
end

function [str] = mesh_face_shading_list
str = verbatim;
%{
                MESH_FACE_SHADING_LIST {
                    %s
                }
%}

function [str] = mesh_face_diffuse_color_list
str = verbatim;
%{
                MESH_FACE_DIFFUSE_COLOR_LIST {
                    %s
                }
%}

function [str] = model_diffuse_color_list
str = verbatim;
%{
                MODEL_DIFFUSE_COLOR_LIST {
                    %s
                }
%}

function [str] = shading_description_str
str = verbatim;
%{
                    SHADING_DESCRIPTION %d {
                         TEXTURE_LAYER_COUNT 0
                         SHADER_ID %d
                    }
%}
