function [vertices, faces, facevertexcdata, renderer, tags] = u3d_pre_patch(ax)
%U3D_PRE_PATCH    Preprocess surface output to u3d.
%
% usage 
%   [vertices, faces, facevertexcdata] = U3D_PRE_PATCH
%   [vertices, faces, facevertexcdata] = U3D_PRE_PATCH(h)
%
% optional input
%   ax = axes object handle
%
% output
%   vertices = row vectors of point positions, as row cell array
%              for multiple surfaces
%            = {1 x #surfaces}
%            = {[#vertices x 3], ... }
%   faces = for each surface triangle face, indices of its 3 vertices,
%           these indices refer to the columns of matrix vertices,
%           as row cell array for multiple surfaces
%         = {1 x #surfaces}
%         = {[#faces x 3], ... }
%   facevertexcdata = RGB color information at each vertex,
%                     as row cell array for multiple surfaces
%                   = {1 x #surfaces}
%                   = {[#vertices x 3], ... }
%
% See also fig2idtf, u3d_pre_line, u3d_pre_surface,
%          u3d_pre_quivergroup, u3d_pre_contourgroup.
%
% File:      u3d_pre_patch.m
% Based on:  u3d_pre by Sven Koerner, koerner(underline)sven(add)gmx.de
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.07.16 - 
% Language:  MATLAB R2012a
% Purpose:   preprocess patch children of axes for u3d export
%
% License to use and modify this code is granted freely to all interested,
% as long as the original author is referenced and attributed as such.
% The original author maintains the right to be solely associated with this work.

%todo
%   if FaceColor == 'none' draw only the edges (i.e., wireframe view)

%% input
if nargin < 1
    sh = findobj('flat', 'type', 'patch');
else
    objs = get(ax, 'Children');
    sh = findobj(objs, 'flat', 'type', 'patch');
end

% Remove dummpy patch
isDummyPatch = arrayfun(@(p) isequal(p.Faces, [1 2 3]) && ...
                             isequal(p.Vertices, [0 1;1 1;0 0]) || ...
                             size(p.Faces, 2) ~= 3, sh);
sh = sh(~isDummyPatch);

if isempty(sh)
    disp('No patch objects found.');
    vertices            = [];
    faces               = [];
    facevertexcdata     = [];
    renderer = [];
    return
end

%% process each patch
N = size(sh, 1); % number of patches
vertices = cell(1, N);
faces = cell(1, N);
facevertexcdata = cell(1, N);
renderer = cell(1, N);
for i=1:N
    fprintf(['\nPreprocessing patch No.', num2str(i), '\n'] );
    h = sh(i, 1);
    
    [v, f, fvx, r] = single_patch_preprocessor(h);
    disp(['Face Vertex Size of patch: ', num2str(size(fvx))]);
    vertices{1, i} = v;
    faces{1, i} = f;
    facevertexcdata{1, i} = fvx;
    renderer{1, i} = r;
end

tags = {sh.Tag};


function [vertices, faces, facevertexcdata, renderer] = single_patch_preprocessor(h)
%% shading -> renderer in adobe reader
edgecolor = get(h, 'EdgeColor');
if strcmp(edgecolor, 'none')
    renderer = 'Solid';
else
    renderer = 'SolidWireframe';
end

%% get defined data-points
faces = get(h, 'Faces');
vertices = get(h, 'Vertices');

facevertexcdata = get(h, 'FaceVertexCData');

% triangulation needed ?
if size(faces, 2) > 3
    faces = delaunay_triangulate_patch_face(faces, vertices);
end

%% remove nan vertices and faces
% remove faces using at least one vertex with some nan coordinate
nan_vertices = any(isnan(vertices), 2);
nan_faces = nan_vertices(faces);
nan_faces = any(nan_faces, 2);
nan_faces = ~nan_faces;
faces = faces(nan_faces, :);

% vertices with nan are not used anymore
% just make them contain numbers
% DO NOT REMOVE them! This would destroy face indexing
vertices(isnan(vertices) ) = 0;
vertices = vertices.';

%% fix face color if needed (after faces have been removed)
facecolor = get(h, 'FaceColor');
if isempty(facevertexcdata) && ~ischar(facecolor)
    disp('Patch: Fixing face color')
    nfaces = size(faces, 1);
    facevertexcdata = repmat(facecolor, nfaces, 1);
end

%% CData
[~, m] = size(facevertexcdata);

% true color ?
if m == 3
    disp('FaceVertexCData is already True Color.')
    return
end

%% continue only if indexed colors need be replaced with RGB from colormap
% not indexed colors ?
if m ~= 1
    error('u3dpatch:FaceVertexCData', 'size(FaceVertexCData, 2) \notin {1, 3}')
end

ax = get(h, 'Parent');

% single row ? - fix to avoid errors later when exporting
% if size(facevertexcdata, 1) 
%     disp('Patch: Single row FaceVerteXCData, replicating for all faces.')
%     nfaces = size(faces, 1);
%     facevertexcdata = repmat(facevertexcdata, nfaces, 1);
% end

facevertexcdata = scaled_ind2rgb(facevertexcdata, ax);


function [realcolor] = scaled_ind2rgb(cdata, ax)
cdata = double(cdata);

cmap = colormap(ax);


realcolor=squeeze(ind2rgb(round(cdata),cmap));


function [f] = delaunay_triangulate_patch_face(faces, vertices)
warning('patch:faces', 'Patch faces are not triangulated.')
disp('      Using Delanay triangulation.')

% for each face
n = size(faces, 1);
f = cell(n, 1);
for i=1:n
    curface = faces(i, :);

    % clear NaN due to face with smaller number of vertices
    curface = curface(~isnan(curface) );

    cur_pntnum = size(curface, 2);

    % triangle after NaN removal ?
    if cur_pntnum == 3
        f(i, 1) = curface;
        continue
    end

    % project on polygon's plane and Delaunay triangulate
    xabs = vertices(curface, :).';
    xref = xabs(:, 1);
    xrel = bsxfun(@minus, xabs, xref);
    base = xrel(:, 2:3); % 1 is the zero vector
    orth_base = orth(base);
    A = orth_base.';
    xproj = A *xrel;
    curtrifaces = delaunay(xproj.');
    
    f{i, 1} = curface(curtrifaces);
end
f = cell2mat(f);
