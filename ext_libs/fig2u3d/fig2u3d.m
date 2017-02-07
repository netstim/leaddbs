function [] = fig2u3d(ax, fname, imgtype, addaxes, varargin)
%FIG2U3D   Convert figure to U3D file.
%   FIG2U3D saves the current axes as a U3D file for inclusion as an
%   interactive 3-dimensional figure within a PDF. Either LaTeX or Adobe
%   Acrobat can be used to embed the U3D file in the PDF.
%
%   A VWS file is also created, which contains the current camera view of
%   the axes saved. This file can be used to set the figure's default view
%   in the PDF to be the same with the open figure window in MATLAB.
%
%   The media9 LaTeX package can import U3D files with their associated VWS
%   files in a PDF document. It can be found here:
%       http://www.ctan.org/tex-archive/macros/latex/contrib/media9
%
%   For PDF readers which do not render 3D figures, it is possible to
%   include an alternative 2D image as a substitute to the 3D object.
%   For conveniency, the script saves a 2D image together with U3D file.
%   File type and other options for exporting this 2D image can be
%   specified as additional arguments.
%
%   Either export_fig (if available) or the standard print function are
%   used for 2D export. Arguments defining 2D export options need conform
%   to the export function used.
%
%   The IDTF2U3D executables are needed from:
%       http://sourceforge.net/projects/u3d/
%   To obtain them download:
%       http://www.mathworks.com/matlabcentral/fileexchange/25383-matlab-mesh-to-pdf-with-3d-interactive-object
%   and place the "bin" directory in the "idtf2u3d" directory of the
%   fig2u3d distribution.
%
% usage
%   FIG2U3D
%   FIG2U3D(ax, fname, imgtype, addaxes, varargin)
%
% optional input
%   ax = axes object handle to export (default = gca)
%   fname = file name string, e.g. 'myfigure' (default = 'surface')
%           Note: Filename extensions are appended automatically.
%   imgtype = save axes also as an image, to be used as a substitute in
%             PDF readers which cannot render the interactive 3D figure.
%           = 'filetype' | 'none'
%             where:
%                   filetype = any file format supported by the export_fig
%                              or print function (whichever is used)
%                   'none' = in this case a 2D image is not saved
%   addaxes = show axes in u3d file (default = 0)
%           = 0 | 1 (do not show/ show, respectively)
%             (depending on which one is available and used).
%   varargin = additional options to be passed to print or export_fig
%              functions, depending on which one is used for saving the 
%              2D image for PDF viewers without 3D support. These must
%              conform to the inputs accepted by the export function used.
%
% output
%   This M-function does not return any data.
%   It saves a U3D file containing the axes object with handle ax.
%
% examples
%   peaks
%   FIG2U3D(gca, 'peaks', '-dpdf') % if using print
%   FIG2U3D(gca, 'peaks', '-dpng', 0, '-r300') % if using print
%
%   FIG2U3D(gca, 'peaks', '-pdf') % if using export_fig
%   FIG2U3D(gca, 'peaks', '-png') % if using export_fig
%   
% remark
%   When the axes are unequal, the exporter maintains the aspect ratio by
%   changing the exported coordinates to achieve the same view as in the
%   MATLAB plot window. To avoid this, edit the "fix_daspect = 1;" line in
%   fig2u3d to "fix_daspect = 0;".
%
% dependency
%   IDTF2U3D executables, download:
%       http://www.mathworks.com/matlabcentral/fileexchange/25383-matlab-mesh-to-pdf-with-3d-interactive-object
%   and place the "bin" directory within the "idtf2u3d" directory of the
%   fig2u3d distribution.
%
% optional dependency
%   export_fig, for saving an accompanying 2D image to substitute the 3D
%   inetractive figure in PDF readers which do not render it.
%   This can be downloaded from the MATLAB Central File Exchange:
%       http://www.mathworks.com/matlabcentral/fileexchange/23629
%   Otherwise, the print function is used for the 2D image.
%
% reference
%   IDTF (Intermediate Data Text File) Format Description, Version 100,
%   Intel Corporation, 2005, available at:
%       http://u3d.svn.sourceforge.net/viewvc/u3d/releases/Gold12Update/Docs/IntermediateFormat/IDTF%20Format%20Description.pdf
%
% See also EXAMPLE_FIG2U3D, FIG2PDF3D, FIG2IDTF, IDTF2U3D, VIEW2VWS, PRINT,
%          EXPORT_FIG.
%
% File:      fig2u3d.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.10 - 2012.07.26
% Language:  MATLAB R2012a
% Purpose:   convert figure to U3D file format (lines, surfaces, points)
% Copyright: Ioannis Filippidis, 2012-

% depends
%   u3d_pre_surface, u3d_pre_line, u3d_pre_quivergroup, u3d_pre_contourgroup
%   fig2idtf, idtf2u3d, view2vws
%   axes_extremal_xyz, takehold, restorehold, quivermd, check_file_extension

% optional depends
%   export_fig

% known limitations
%   surface texturemapping not exported
%   filled contourplots not exported

% todo
%   use viewmtx

%% input
if nargin < 1
    ax = gca;
else
    % not a handle ?
    if ~ishandle(ax)
        warning('fig2u3d:handle',...
                'AX argument is not a handle. Exporting current axes.')
        ax = gca;
    end
    
    % not axes handle ?
    handle_type = get(ax, 'Type');
    if strcmp(handle_type, 'figure')
        warning('fig2u3d:handle',...
               ['AX argument is a figure handle.',...
                'Exporting the first axes child of this figure.'] )
        ax = get(ax, 'Children');
    end
end

% no axes object ?
if isempty(ax)
    warning('fig2u3d:emptyaxes',...
            'AX axes handle is empty. Noting to export.')
    return
end

% filename provided ?
if nargin < 2
    fname = 'surface';
elseif ~ischar(fname) || isempty(fname)
    fname = 'surface';
end

% save image substitute ?
saveimg = 1;
if nargin < 3
    imgtype = '-png';
elseif strcmp(imgtype, 'none') || isempty(imgtype)
    saveimg = 0;
elseif ~ischar(imgtype)
    imgtype = '-png';
else
    imgtype = '-pdf';
end

% show axes ?
if nargin < 4
    addaxes = 0;
elseif isempty(addaxes)
    addaxes = 0;
elseif (addaxes ~= 1) && (addaxes ~= 0)
    error('fig2u3d:addaxes', 'Argument "addaxes" can be 0 or 1.')
end

plot_axes(ax, addaxes)
delete_idtf = 1; % delete intermediate IDTF file

% nothing in the plot ?
obj = get(ax, 'Children');
if isempty(obj)
    warning('axes:empty', 'Nothing in the plot. U3D file not exported.')
    return
end

% are the axes equal ?
dar = daspect(ax);
dar = dar /dar(1);
if any(dar ~= 1)
    msg = 'Axis aspect ratio unequal. Exported U3D will look different.';
    warning('axes:aspect', msg)
end

fix_daspect = 1; % when axes are unequal, change aspect ratio of exported
                 % graphics, so that the exported figure has the same view

%% convert graphics objects to meshes, line_sets and point_sets
[surf_vertices, surf_faces, surf_facevertexcdata, surf_renderers] = u3d_pre_surface(ax);
[patch_vertices, patch_faces, patch_facevertexcdata, patch_renderers] = u3d_pre_patch(ax);
[line_vertices, line_edges, line_colors,...
                line_points, line_point_colors] = u3d_pre_line(ax);
[quiver_vertices, quiver_edges, quiver_colors] = u3d_pre_quivergroup(ax);
[contour_vertices, contour_edges, contour_colors] = u3d_pre_contourgroup(ax);

%% group meshes, line_sets and point_sets
% aggregate meshes
mesh_vertices = [surf_vertices, patch_vertices];
mesh_faces = [surf_faces, patch_faces];
mesh_colors = [surf_facevertexcdata, patch_facevertexcdata];

% aggregate lines
line_vertices = [line_vertices, quiver_vertices, contour_vertices];
line_edges = [line_edges, quiver_edges, contour_edges];
line_colors = [line_colors, quiver_colors, contour_colors];

% aggregate pointsets
pointset_points = line_points;
pointset_colors = line_point_colors;

%% fix axis aspect ratio
if fix_daspect == 1
    dar = daspect(ax);
    
    if ~isequal(dar, ones(1, 3) )
        msg = ['Unequal axes data aspect ratio, to be fixed by ',...
               'scaling the axes. To change this behavior, edit ',...
               '"fix_daspect = 1;" in fig2u3d to "fix_daspect = 0;"'];
        warning('fig2u3d:daspect', msg)
    end
    
    mesh_vertices = axis_rescale(mesh_vertices, dar);
    line_vertices = axis_rescale(line_vertices, dar);
    pointset_points = axis_rescale(pointset_points, dar);
end

%% export

fig2idtf(fname,...
          mesh_vertices, mesh_faces, mesh_colors,...
          line_vertices, line_edges, line_colors,...
          pointset_points, pointset_colors)

idtf2u3d(fname)
rm_idtf(fname, delete_idtf)

part_renderers = [surf_renderers, patch_renderers];
view2vws(ax, fname, part_renderers, fix_daspect)

save_png_substitute(ax, fname, saveimg, imgtype, varargin{:} )

function rm_idtf(fname, yesno)
if yesno == 0
    return
end

% fname extensions ok ?
fname = check_file_extension(fname, '.idtf');
fname  = fullfile(cd, fname);

system(['rm ', fname] )

function plot_axes(ax, addaxes)
if addaxes == 0
    return
end

[xyz_minmax, dim] = axes_extremal_xyz(ax);
xyz_minmax = 1.3 *xyz_minmax;

x0 = [0, 0, 0].';
vx = [xyz_minmax(2), 0, 0].';
vy = [0, xyz_minmax(4), 0].';

if dim == 3
    vz = [0, 0, xyz_minmax(6) ].';
end

held = takehold(ax);
    quivermd(ax, x0, vx, 'k')
    quivermd(ax, x0, vy, 'k')
    quivermd(ax, x0, vz, 'k')
restorehold(ax, held)

function [] = save_png_substitute(ax, fname, saveimg, imgtype, varargin)
% save image substitute ?
if saveimg == 0
    return
end

% export_fig available ?
fighandle = get(ax, 'Parent');
if exist('export_fig', 'file') == 2
    col = get(fighandle, 'Color'); % current fig background color
    
    % set white background for png ?
    if strcmp(imgtype, '-png')
        set(fighandle, 'Color', 'w')
    end
    
    export_fig(ax, imgtype, varargin{:}, fname)
    
    set(fighandle, 'Color', col) % restore color
else
    if strfind(imgtype, 'pdf')
        imgtype = '-dpdf';
    elseif strfind(imgtype, 'png')
        imgtype = '-dpng';
    end
    print(fighandle, imgtype, varargin{:}, fname)
end

function [a] = axis_rescale(a, dar)
if isempty(a)
   return
end

a = cellfun(@(x) cell_divider(x, dar), a, 'UniformOutput', false);

function [x] = cell_divider(x, dar)
x = [x(1, :)/dar(1);
     x(2, :)/dar(2);
     x(3, :)/dar(3) ];
