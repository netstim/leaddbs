function [line_resources] = populate_line_resource_str(line_vertices,...
                                line_lines, line_colors, n_resources)
%
% See also FACE_VERTEX_DATA_EQUALS_NPOINTS, VERBATIM.
%
% File:      populate_line_resource_str.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.11 - 
% Language:  MATLAB R2012a
% Purpose:   line resource strings
% Copyright: Ioannis Filippidis, 2012-

% depends
%   verbatim

n_lines = size(line_vertices, 2);
tmp_line_resources = cell(1, n_lines);
for i=1:n_lines
    cur_line_points = line_vertices{1, i};
    cur_line_edges = line_lines{1, i};
    cur_line_color = line_colors{1, i};
    
    cur_line_resource = single_line_resource_str(cur_line_points,...
                         cur_line_edges, i, cur_line_color, n_resources);
    
    tmp_line_resources{1, i} = ['\n\n', cur_line_resource];
end
line_resources = [tmp_line_resources{:} ];
%disp(line_resources)

function [line_resource_str] = single_line_resource_str(vertices, edges,...
                                    line_number, line_color, n_resources)
npnt = size(vertices, 2);
nlines = size(edges, 2);

% at least one line ?
% (needs 2 vertices)
if nlines < 1
    msg = ['No edges for this Lineset (id = ', num2str(line_number), ' )'];
    error('idtf:line', msg)
end

line_position_list = edges;

str_line_position_list = sprintf('    %d %d\n', line_position_list);
str_line_normal_list = sprintf('    %d %d\n', [1:nlines; 1:nlines]-1);
str_line_shading_list = sprintf('    %d\n', zeros(nlines, 1) );
str_line_diffuse_color_list = sprintf('%d\n', zeros(nlines, 1) );

str_model_position_list = sprintf('    %1.6f %1.6f %1.6f\n', vertices);

model_normal_count = nlines;
str_model_normal_list = sprintf('    %f %f %f\n', zeros(3, nlines) );

resource_number = n_resources +line_number -1;

% fliplr is used because of a bug in IDTFConverter.exe
% which recognizes BGR instead of RGB
model_diffuse_color_list = sprintf('    %f %f %f\n', fliplr(line_color) );

line_resource_str = sprintf(resource_line_str,...
    resource_number, line_number,...
    nlines, npnt, model_normal_count,...
    str_line_position_list, str_line_normal_list, str_line_shading_list,...
    str_line_diffuse_color_list,...
    str_model_position_list, str_model_normal_list, model_diffuse_color_list);

function [str] = resource_line_str
% resource_number, line_number, line_count, model_position_count,
% model_normal_count,
% line_position_list, line_normal_list, line_shading_list,
% model_position_list, model_normal_list, model_diffuse_color_list
str = verbatim;
%{

RESOURCE %d {
    RESOURCE_NAME "MyLine%d"
    MODEL_TYPE "LINE_SET"
    LINE_SET {
        LINE_COUNT %d
        MODEL_POSITION_COUNT %d
        MODEL_NORMAL_COUNT %d
        MODEL_DIFFUSE_COLOR_COUNT 1
        MODEL_SPECULAR_COLOR_COUNT 0
        MODEL_TEXTURE_COORD_COUNT 0
        MODEL_SHADING_COUNT 1
        MODEL_SHADING_DESCRIPTION_LIST {
            SHADING_DESCRIPTION 0 {
                TEXTURE_LAYER_COUNT 0
                SHADER_ID 0
            }
        }
        LINE_POSITION_LIST {
            %s
        }
        LINE_NORMAL_LIST {
            %s
        }
        LINE_SHADING_LIST {
            %s
        }
        LINE_DIFFUSE_COLOR_LIST {
            %s
        }
        
        MODEL_POSITION_LIST {
            %s
        }
        MODEL_NORMAL_LIST {
            %s
        }
        MODEL_DIFFUSE_COLOR_LIST {
            %s
        }
    }
}

%}
