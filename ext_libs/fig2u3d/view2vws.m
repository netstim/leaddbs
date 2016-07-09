function [] = view2vws(ax, filename, part_renderers, fix_daspect)
%VIE2VWS    Saves current view in a views file for LaTeX media9 package.
%
% See also FIG2U3D.
%
% File:      view2vws.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.10 - 2012.07.16
% Language:  MATLAB R2012a
% Purpose:   export view for LaTeX media9 package
% Copyright: Ioannis Filippidis, 2012-

% depends
%   axes_extremal_xyz, verbatim, check_file_extension

%% input
if nargin < 1
    ax = gca;
end

% filename and its extension ?
if nargin < 2
    filename = 'matlab.vws';
else
    filename = check_file_extension(filename, '.vws');
end
matrix_mode = 0; % matrix mode not implemented yet
viewname = 'MATLABfig';

% part renderers
if nargin < 3
    part_renderers = [];
end

if nargin < 4
    fix_daspect = 0;
end

%% get view
camera_position = get(ax, 'CameraPosition'); % = campos

center_of_orbit = get(ax, 'CameraTarget'); % 3Dcoo
center_of_orbit_2_camera_vector = camera_position -center_of_orbit; % 3Dc2c
radius_of_orbit = norm(center_of_orbit_2_camera_vector); % 3Droo
camera_roll = 0; % 3Droll

% [3x4] homogenous transformation (sub)matrix
curview = view(ax);
T = curview(1:3, :); % reduced, 12-element matrix for Adobe PDF internals
camera_2_world_transformation_matrix = T(:).'; % c2w

camera_aperture_angle = get(ax, 'CameraViewAngle'); % 3Daac (degrees) = camva

%% data aspect ratio fixing
if fix_daspect == 1
    dar = daspect(ax);
    %length_ratio = 3^0.5 /norm(dar); % max (new length) /(old length)
    
    center_of_orbit = center_of_orbit ./dar;
    center_of_orbit_2_camera_vector = center_of_orbit_2_camera_vector ./dar;
    %radius_of_orbit = radius_of_orbit;
    %orthographic_scaling_factor = orthographic_scaling_factor;
end

projection = get(ax, 'Projection');
if strcmp(projection, 'orthographic')
    xyz_minmax = axes_extremal_xyz(ax);
    
    d1 = xyz_minmax(2) -xyz_minmax(1);
    d2 = xyz_minmax(4) -xyz_minmax(3);
    d3 = xyz_minmax(6) -xyz_minmax(5);
    
    d123 = [d1, d2, d3];
    
    if fix_daspect == 1
        d123 = d123 ./dar;
    end
    
    d = norm(d123); % minimal bounding sphere radius
    
    orthographic_scaling_factor = 1 /d; % 3Dortho (scaling factor)
    radius_of_orbit = d /2; % 3Droo
    
    msg = ['Orthographic projection: Camera distance to target is different',...
           ' from figure. Computed based on the scene bounding sphere.'];
    warning('view:export', msg)
else
    orthographic_scaling_factor = 0;
end

%% views
fileinfo = verbatim;
%{
%% File type:    VWS
%% Version:      2012-05-11
%% Created by:   fig2u3d, MATLAB figure export tool by Ioannis Filippidis
%% Help:         See LaTeX package media9 documentation by Alexander Grahn
%% Description:  View settings file for u3d file to be included in pdf
\n
%}
viewname = ['VIEW=', viewname, '\n'];
closing = 'END\n';

if orthographic_scaling_factor ~= 0
    proj = ['    ORTHO=', num2str(orthographic_scaling_factor), '\n'];
else
    proj = ['    AAC=', num2str(camera_aperture_angle), '\n'];
end

if matrix_mode == 1
    c2w = ['    C2W=', num2str(camera_2_world_transformation_matrix), '\n'];
    
    views = [viewname, c2w, proj];
else
    coo = ['    COO=', num2str(center_of_orbit), '\n'];
    c2c = ['    C2C=', num2str(center_of_orbit_2_camera_vector), '\n'];
    roo = ['    ROO=', num2str(radius_of_orbit), '\n'];
    
    roll = ['    ROLL=', num2str(camera_roll), '\n'];
    
    views = [viewname, coo, c2c, roo, roll, proj];
end

str = [fileinfo, views];

%% part renderers
if ~isempty(part_renderers)
    part_strs = add_part_renderers(part_renderers);
else
    part_strs = {''};
end

str = [str, part_strs{:} ];
str = [str, closing];

%% export to file
fid = fopen(filename, 'w');
fprintf(fid, str);
fclose(fid);

%% screen msg
s = sprintf(str);
disp('Exported axes view in .vws file for LaTeX media9 package is:')
disp(s)

function [part_strs] = add_part_renderers(part_renderers)
n = size(part_renderers, 2);
part_strs = cell(1, n);
for i=1:n
    s = part_str;
    
    partname = ['Mesh', num2str(i) ];
    part_rendermode = part_renderers{1, i};
    
    s = sprintf(s, partname, part_rendermode);
    part_strs{1, i} = s;
end

function [str] = part_str
str = verbatim;
%{

    PART=%s
        RENDERMODE=%s
    END

%}
