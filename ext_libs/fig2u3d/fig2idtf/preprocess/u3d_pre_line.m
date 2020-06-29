function [vertices, edges, colors, points, point_colors] = u3d_pre_line(ax)
%U3D_PRE_LINE   Preprocess line output to u3d.
%
% usage
%   points = U3D_PRE_LINE
%   points = U3D_PRE_LINE(ax)
%
% optional input
%   ax = axes object handle
%
% output
%   points = position vectors as columns of matrix, as row cell array,
%            for multiple lines
%          = {1 x #lines}
%          = {[3 x #points], ... }
%
% See also fig2idtf, u3d_pre_surface, u3d_pre_patch, u3d_pre_quivergroup,
%          u3d_pre_contourgroup.
%
% File:      u3d_pre_line.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.10 - 
% Language:  MATLAB R2012a
% Purpose:   preprocess lineseries children of axes for u3d export
% Copyright: Ioannis Filippidis, 2012-

% depends
%   create_marker_lines, cut_line_to_pieces, get_line_xyz
%   arclength, meshgrid2vec

%% input
if nargin < 1
    sh = findobj('type', 'line');
else
    objs = get(ax, 'Children');
    % using full depth may be advantageous, but nans should then be treated instead
    sh = findobj(objs, 'flat', 'type', 'line');
end

if isempty(sh)
    vertices = [];
    edges = [];
    colors = [];
    
    points = [];
    point_colors = [];
    disp('No lines found.');
    return
end

%% process each line
N = size(sh, 1); % number of lines
vertices = cell(0); % init is pointless when chopping is implemented
edges = cell(0);
colors = cell(0);
k = 0;
points = cell(0);
point_colors = cell(0);
m = 0;
for i=1:N
    disp(['     Preprocessing line No.', num2str(i) ] );
    h = sh(i, 1);
    
    % linestyle represented by 
    [curvertices, curedges, curcolor] = line_prep(h);
    
    % add linesets due to linestyle
    n = size(curvertices, 2);
    if n ~= 0
        I = k +(1:n);

        vertices(1, I) = curvertices;
        edges(1, I) = curedges;
        colors(1, I) = curcolor;
    
        k = k +n;
    end
    disp(['         used: ', num2str(n), ' LineSets for LineStyle.'] )
    
    [curvertices, curedges, curcolor, curpoints, curpoint_colors] = marker_prep(h);
    
    % add linesets due to markers
    n = size(curvertices, 2);
    if n ~= 0
        I = k +(1:n);

        vertices(1, I) = curvertices;
        edges(1, I) = curedges;    
        colors(1, I) = curcolor;

        k = k +n;
    end
    disp(['         used: ', num2str(n), ' LineSets for MarkerStyle.'] )
    
    % add pointsets due to linesets and markers
    n = size(curpoints, 2);
    if n ~= 0
        I = m +(1:n);
    
        points(1, I) = curpoints;
        point_colors(1, I) = curpoint_colors;

        m = m +n;
    end
    disp(['         used: ', num2str(n), ' PointSets for MarkerStyle.'] )
    
    disp(['         cumulative at: ', num2str(k), ' LineSets.'] )
    disp(['         cumulative at: ', num2str(m), ' PointSets.'] )
    
    %vertices{1, i} = curvertices;
    %edges{1, i} = curedges;
end

function [vertices, edges, line_colors, points, point_color] = marker_prep(h)
% init linesets
vertices = {};
edges = {};
line_colors = {};

% init pointsets
points1 = {};
point_color1 = {};

points2 = {};
point_color2 = {};

%% interpolate dots ?
linestyle = get(h, 'LineStyle');
if strcmp(linestyle, ':')
    % (this is linestyle producing points)
    % interpolate
    p = get_line_xyz(h);

    points1 = dots_linestyle(p);
    
    n = size(points1, 2);
    point_color1 = {get(h, 'Color') };
    point_color1 = repmat(point_color1, 1, n);
end

%% Marker Style ?
markerstyle = get(h, 'Marker');
switch markerstyle
    case 'none'
        disp('         This line has: Marker = ''none''.')
    case '.'
        % . = point
        % marker producing pointset
        points2 = get_line_xyz(h);
        
        % reduce problems due to compression by IDTFConverter
        piece_size = 10;
        points2 = cut_line_to_pieces(points2, piece_size);
        
        n = size(points2, 2);
        point_color2 = {get(h, 'Color') };
        point_color2 = repmat(point_color2, 1, n);
    otherwise
        [vertices, edges, line_colors] = create_marker_lines(h, markerstyle);
end

points = [points1, points2];
point_color = [point_color1, point_color2];

function [points] = dots_linestyle(p)
[L, dL] = arclength(p);
d = L /100; % distance between dots

dL = dL.'; % [1 x #segments]
n = size(dL, 2); % number of segments

points = cell(1, n);
for i=1:n
    p1 = p(:, i);
    p2 = p(:, i+1);
    
    curdL = dL(1, i);
    
    m = floor(curdL /d) +1; % # of dots between two consecutive line points
    dots = linspacendim(p1, p2, m);
    
    points{1, i} = dots;
end

function [vertices, edges, line_colors] = line_prep(h)
% no line visible ?
linestyle = get(h, 'LineStyle');
disp(['         This line has: LineStyle = ''', linestyle, '''.'] )
switch linestyle
    case {'none', '.', ':'}
        disp('      No line for this lineseries exported to u3d.')
        vertices = {};
        edges = {};
        line_colors = {};
        return
    case '-'
        % solid or dashdot
        % (the dot part of dashdot is handled by marker_prep)
        dashratio = 0;
        [vertices, edges, line_colors] = solid_line_prep(h, dashratio);
    case '--'
        % dashed
        dashratio = 0.5;
        [vertices, edges, line_colors] = solid_line_prep(h, dashratio);
    otherwise
        error(['Unsupported LineStyle:', linestyle] )
end

function [vertices, edges, line_colors] = solid_line_prep(h, dashratio)
lv = get_line_xyz(h); % line evrtices
line_color = get(h, 'Color');

% single point only ?
npoints = size(lv, 2);
if npoints == 1
    vertices = {};
    edges = {};
    line_colors = {};
    return
end

%edges = [1:(npnt-1); 2:npnt]-1;

%% temporary solution to reduces compression problems after idtf2u3d.exe
piece_size = 10;
[vertices, edges] = line_pieces(lv, piece_size, dashratio);
n = size(vertices, 2);
line_colors = repmat({line_color}, 1, n);
