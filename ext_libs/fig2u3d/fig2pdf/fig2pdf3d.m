function [] = fig2pdf3d(ax, filename, media9_or_movie15, pdforxelatex)
%FIG2PDF3D  Convert axes to PDF with embedded interactive 3D image.
%
% usage
%   FIG2PDF3D(ax, filename, media9_or_movie15, pdforxelatex)
%
% input
%   ax = axes object handle
%   filename = file name string (default = 'surface')
%   media9_or_movie15 = use media9 or movie15 LaTeX packages
%                     = 'media9' | 'movie15' (default = 'media9')
%   pdforxelatex = use pdf or xelatex
%                = 'pdflatex' | 'xelatex'
%
% output
%   Saves the axes object as a PDF with interactive 3D content.
%
% remark
%   The movie15 LaTeX package implements only perspective projection
%
% dependency
%   pdfLaTeX or XeLaTeX
%   media9 or movie15 LaTeX packages
%
% See also FIG2LATEX, LATEX2PDF3D, FIG2U3D, FIG2IDTF, IDTF2U3D.
%
% File:      fig2pdf3d.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.22
% Language:  MATLAB R2012a
% Purpose:   convert MATLAB figure to interactive 3D image embedded in PDF
%
% acknowledgment
%   Based on demo_mesh2pdf by Alexandre Gramfort.
%   This can be found on the MATLAb Central File Exchange:
%       http://www.mathworks.com/matlabcentral/fileexchange/25383-matlab-mesh-to-pdf-with-3d-interactive-object
%   and is covered by the BSD License.

% depends
%   fig2latex, latex2pdf3d, {pdflatex | xelatex}, {media9 | movie15}

%% input
if nargin < 1
    ax = gca;
end

% Set the output Latex filename
if nargin < 2
    filename = 'surface';
end

% which package to use ?
if nargin < 3
    media9_or_movie15 = 'media9';
end

% which LaTeX compiler ?
if nargin < 4
    pdforxelatex = 'xelatex';
end

%% Generate Latex file
fig2latex(ax, filename, media9_or_movie15, pdforxelatex);
latex2pdf3d(filename, pdforxelatex)
rm_aux_files(filename)

function [] = rm_aux_files(fname)
fname = clear_file_extension(fname, '.tex');

delete([fname, '.png'] )
delete([fname, '.tex'] )
delete([fname, '_small.tex'] )
delete([fname, '.u3d'] )
delete([fname, '.vws'] )
delete([fname, '.aux'] )
delete([fname, '.log'] )
