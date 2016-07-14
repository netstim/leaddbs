function [pointset_resources] = populate_point_resource_str(points, colors, n_resources)
%
% See also VERBATIM.
%
% File:      populate_point_resource_str.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.14 - 
% Language:  MATLAB R2012a
% Purpose:   point resource strings
% Copyright: Ioannis Filippidis, 2012-

% depends
%   verbatim

n_pointsets = size(points, 2);
pnt_rc = cell(1, n_pointsets);
for i=1:n_pointsets
    cur_points = points{1, i};
    cur_color = colors{1, i};
    
    cur_resource = single_pointset_resource_str(cur_points, cur_color,...
                                                i, n_resources);
    
    pnt_rc{1, i} = ['\n\n', cur_resource];
end
pointset_resources = [pnt_rc{:} ];
%disp(pointset_resources)

function [line_resource_str] = single_pointset_resource_str(points, point_color,...
                                    pointset_id, n_resources)
str = pointset_resource_str;
npnt = size(points, 2);

resource_id = n_resources +pointset_id -1;

point_position_list = 0:(npnt-1);
point_normal_list = 0:(npnt-1);

point_position_list_str = sprintf('%d\n', point_position_list);
point_normal_list_str = sprintf('%d\n', point_normal_list);
point_shading_list_str = sprintf('%d\n', zeros(npnt, 1) );
point_diffuse_color_list_str = sprintf('%d\n', zeros(npnt, 1) );

model_position_list_str = sprintf('%1.6f %1.6f %1.6f\n', points);

model_normal_count = npnt;
model_normal_list_str = sprintf('%f %f %f\n', zeros(3, npnt) );

% fliplr is used because of a bug in IDTFConverter.exe
% which recognizes BGR instead of RGB
model_diffuse_color_list_str = sprintf('%f %f %f\n', fliplr(point_color) );

line_resource_str = sprintf(str,...
    resource_id, pointset_id, npnt,...
    npnt, model_normal_count,...
    point_position_list_str, point_normal_list_str,...
    point_shading_list_str, point_diffuse_color_list_str,...
    model_position_list_str, model_normal_list_str, model_diffuse_color_list_str);

function [str] = pointset_resource_str
% resource_id, pointset_id, point_count,
%
% model_position_count, model_normal_count,
%
% point_position_list, point_normal_list,
% point_shading_list, point_diffuse_color_list,
%
% model_position_list, model_normal_list, model_diffuse_color_list
str = verbatim;
%{

RESOURCE %d {
    RESOURCE_NAME "MyPoints%d"
    MODEL_TYPE "POINT_SET"
    POINT_SET {
        POINT_COUNT %d
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
        POINT_POSITION_LIST {
            %s
        }
        POINT_NORMAL_LIST {
            %s
        }
        POINT_SHADING_LIST {
            %s
        }
        POINT_DIFFUSE_COLOR_LIST {
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
