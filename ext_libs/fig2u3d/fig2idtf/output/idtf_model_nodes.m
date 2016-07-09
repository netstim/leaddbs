function [nodes] = idtf_model_nodes(n_meshes, n_lines, n_pointsets)
% File:      idtf_model_nodes.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.10 - 2012.06.24
% Language:  MATLAB R2012a
% Purpose:   create IDTF model nodes
% Copyright: Ioannis Filippidis, 2012-

% depends
%   verbatim

%% mesh nodes
mesh_node = mesh_node_str;
temp_nodes = cell(1, n_meshes);
for i=1:n_meshes
    cur_mesh_node = sprintf(mesh_node, i, i);
    temp_nodes{1, i} = cur_mesh_node;
end
mesh_nodes = [temp_nodes{:} ];
%disp(mesh_nodes)

%% line nodes
line_node = line_node_str;
temp_nodes = cell(1, n_lines);
for i=1:n_lines
    cur_line_node = sprintf(line_node, i, i);
    temp_nodes{1, i} = cur_line_node;
end
line_nodes = [temp_nodes{:} ];
%disp(line_nodes)

%% pointset nodes
point_node = point_node_str;
temp_nodes = cell(1, n_pointsets);
for i=1:n_pointsets
    cur_point_node = sprintf(point_node, i, i);
    temp_nodes{1, i} = cur_point_node;
end
point_nodes = [temp_nodes{:} ];
%disp(point_nodes)

%% output
nodes = [mesh_nodes, line_nodes, point_nodes];

function [str] = mesh_node_str
%mesh_id_number, mesh_id_number
str = verbatim;
%{

NODE "MODEL" {
     NODE_NAME "Mesh%d"
     PARENT_LIST {
          PARENT_COUNT 1
          PARENT 0 {
               PARENT_NAME "<NULL>"
               PARENT_TM {
                    1.000000 0.000000 0.000000 0.000000
                    0.000000 1.000000 0.000000 0.000000
                    0.000000 0.000000 1.000000 0.000000
                    0.000000 0.000000 0.000000 1.000000
               }
          }
     }
     RESOURCE_NAME "MyMesh%d"
     MODEL_VISIBILITY "BOTH"
}

%}

function [str] = line_node_str
% lineset_id_number, lineset_id_number
str = verbatim;
%{

NODE "MODEL" {
     NODE_NAME "Line%d"
     PARENT_LIST {
          PARENT_COUNT 1
          PARENT 0 {
               PARENT_NAME "<NULL>"
               PARENT_TM {
                    1.000000 0.000000 0.000000 0.000000
                    0.000000 1.000000 0.000000 0.000000
                    0.000000 0.000000 1.000000 0.000000
                    0.000000 0.000000 0.000000 1.000000
               }
          }
     }
     RESOURCE_NAME "MyLine%d"
}

%}

function [str] = point_node_str
% pointset_id_number, pointset_id_number
str = verbatim;
%{

NODE "MODEL" {
     NODE_NAME "Points%d"
     PARENT_LIST {
          PARENT_COUNT 1
          PARENT 0 {
               PARENT_NAME "<NULL>"
               PARENT_TM {
                    1.000000 0.000000 0.000000 0.000000
                    0.000000 1.000000 0.000000 0.000000
                    0.000000 0.000000 1.000000 0.000000
                    0.000000 0.000000 0.000000 1.000000
               }
          }
     }
     RESOURCE_NAME "MyPoints%d"
}

%}
