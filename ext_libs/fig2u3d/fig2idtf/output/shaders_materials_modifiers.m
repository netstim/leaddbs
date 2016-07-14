function [shaders, materials, modifiers] = shaders_materials_modifiers(v, f, c)
% File:      shaders_materials_modifiers.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.24
% Language:  MATLAB R2012a
% Purpose:   mesh shaders, materials and modifiers
% Copyright: Ioannis Filippidis, 2012-

% depends
%   verbatim

n_meshes = size(v, 2);
n_shaders = 0;
shader_resources = cell(1, n_meshes);
n_materials = 0;
material_resources = cell(1, n_meshes);
modifiers = cell(1, n_meshes);
for i=1:n_meshes
    vertices = v{1, i}.';
    faces = f{1, i};
    face_vertex_data = c{1, i};
    
    npoints = size(vertices, 1);
    nfaces = size(faces, 1);
    nfvd = size(face_vertex_data, 1);
    
    face_vertex_data_unique = unique(face_vertex_data, 'rows');
    nface_vertex_data_unique = size(face_vertex_data_unique, 1);
    
    %% shader and material
    if nfvd == npoints
        use_vertex_color = 1;
        shidx = n_shaders; %0
        
        str2 = sprintf(material_resource_str(use_vertex_color), [shidx, shidx].');
        
        cur_n_shaders = 1;
    elseif nfvd == nfaces
        use_vertex_color = 0;
        shidx = (0:(nface_vertex_data_unique -1) ).' +n_shaders;
        
        str2 = sprintf(material_resource_str(use_vertex_color), [shidx, shidx, face_vertex_data_unique].');
        
        cur_n_shaders = nface_vertex_data_unique;
    else
        warning('shader:color', 'use_vertex_color: nfvd \notin {npoints, nfaces}')
    end
    
    shader_resources{1, i} = sprintf(shader_resource_str(use_vertex_color), [shidx, shidx, shidx].');
    material_resources{1, i} = str2;
    
    n_shaders = n_shaders +cur_n_shaders;
    n_materials = n_materials +cur_n_shaders;
    
    %% modifier
    shader_list_id_and_name = [(0:(cur_n_shaders-1) ).', shidx].';
    str = sprintf(shader_list_str, shader_list_id_and_name);
    modifiers{1, i} = sprintf(shading_modifiers_heading, i, cur_n_shaders, str);
end

%% output
shader_resources = [shader_resources{:} ];
shaders = sprintf(shader_resources_header, n_shaders, shader_resources);

material_resources = [material_resources{:} ];
materials = sprintf(material_resources_header, n_materials, material_resources);

modifiers = [modifiers{:} ];

function [str] = shader_resources_header
str = verbatim;
%{

RESOURCE_LIST "SHADER" {
     RESOURCE_COUNT %d
%s
}

%}

function [str] = shader_resource_str(use_vertex_color)
if use_vertex_color
    str = verbatim;
%{
     RESOURCE %d {
          RESOURCE_NAME "Box01%d"
          ATTRIBUTE_USE_VERTEX_COLOR "TRUE"
          SHADER_MATERIAL_NAME "Box01%d"
          SHADER_ACTIVE_TEXTURE_COUNT 0
    }

%}
else
    str = verbatim;
%{
     RESOURCE %d {
          RESOURCE_NAME "Box01%d"
          SHADER_MATERIAL_NAME "Box01%d"
          SHADER_ACTIVE_TEXTURE_COUNT 0
    }

%}
end

function [str] = material_resources_header
str = verbatim;
%{

RESOURCE_LIST "MATERIAL" {
     RESOURCE_COUNT %d
%s
}

%}

function [str] = material_resource_str(use_vertex_color)
% resource_id, box_id
% OR
% resource_id, box_id, material_diffuse_colors(3)

if use_vertex_color
    str = verbatim; %use for alpha: MATERIAL_OPACITY 0.500000
%{
     RESOURCE %d {
          RESOURCE_NAME "Box01%d"
          MATERIAL_AMBIENT 0.0 0.0 0.0
          MATERIAL_DIFFUSE 1.0 1.0 1.0
          MATERIAL_SPECULAR 0.0 0.0 0.0
          MATERIAL_EMISSIVE 1.0 1.0 1.0
          MATERIAL_REFLECTIVITY 0.500000
          MATERIAL_OPACITY 1.0
     }

%}
else
    str = verbatim;
%{
     RESOURCE %d {
          RESOURCE_NAME "Box01%d"
          MATERIAL_AMBIENT 0.0 0.0 0.0
          MATERIAL_DIFFUSE %f %f %f
          MATERIAL_SPECULAR 0.2 0.2 0.2
          MATERIAL_EMISSIVE 0.0 0.0 0.0
          MATERIAL_REFLECTIVITY 0.100000
          MATERIAL_OPACITY 1.000000
     }

%}
end

function [str] = shading_modifiers_heading
str = verbatim;
%{

MODIFIER "SHADING" {
     MODIFIER_NAME "Mesh%d"
     PARAMETERS {
          SHADER_LIST_COUNT %d
          SHADER_LIST_LIST {
%s
          }
     }
}

%}

function [str] = shader_list_str
str = verbatim;
%{
               SHADER_LIST %d {
                    SHADER_COUNT 1
                    SHADER_NAME_LIST {
                         SHADER 0 NAME: "Box01%d"
                    }
               }

%}
