function [] = fig2latex(ax, fname, media9_or_movie15, pdforxelatex)
%FIG2LATEX  Convert axes to U3D file, generating LaTeX code including it.
%
% usage
%   FIG2LATEX(ax, fname, media9_or_movie15, pdforxelatex)
%
% input
%   ax = axes object handle
%   fname = file name
%   media9_or_movie15 = select LaTeX package for including the U3D file
%                     = 'media9' | 'movie15'
%   pdforxelatex = select LaTeX compiler
%                = 'xelatex' | 'pdflatex'
%
% output
%   Saves a U3D file and the LaTeX which includes it.
%
% remark
%   The movie15 LaTeX package implements only perspective projection
%
% See also FIG2PDF3D, U3D_IN_LATEX.
%
% File:      fig2latex.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.22
% Language:  MATLAB R2012a
% Purpose:   convert MATLAB figure to U3D file and LaTeX code including it
%
% acknowledgment
%   Based on mesh_to_latex by Alexandre Gramfort.
%   This can be found on the MATLAb Central File Exchange:
%       http://www.mathworks.com/matlabcentral/fileexchange/25383-matlab-mesh-to-pdf-with-3d-interactive-object
%   and is covered by the BSD License.

% depends
%   fig2u3d, u3d_in_latex, clear_file_extension, check_file_extension

% change to perspective projection for movie15 ?
if strcmp(media9_or_movie15, 'movie15')
    msg = 'The movie15 LaTeX package only uses Perspective projection.';
    warning('movie15:proj', msg)
    disp('      Changing to ''Projection'' = ''perspective'' ');
    curproj = get(ax, 'Projection');
    set(ax, 'Projection', 'perspective')
    fig2u3d(ax, fname)
    set(ax, 'Projection', curproj)
else
    fig2u3d(ax, fname)
end

u3d_in_latex(fname, media9_or_movie15)

fname = clear_file_extension(fname, '.tex');
texfile  = check_file_extension(fname, '.tex');
disp(['Writing: ', texfile] );

[~, fname] = fileparts(fname);

% XeLaTeX or PDFLaTeX ?
if strcmp(pdforxelatex, 'xelatex')
    s = xelatex;
else
    s = '';
end

doc = '\\documentclass[class=minimal,border=0pt]{standalone}\n';
if strcmp(media9_or_movie15, 'media9')
    content = [doc, s, media9(pdforxelatex), latex_content];
elseif strcmp(pdforxelatex, 'xelatex')
    error('xelatex:movie15', 'Movie15 does not work well with XeLaTeX.')
else
    content = [doc, s, movie15, latex_content];
end

fid = fopen(texfile, 'w');
    fprintf(fid, content, fname);
fclose(fid);

function [str] = latex_content
% media9 or movie15 = string,
% packages needed by xelatex = string
str = verbatim;
%{
\\usepackage{graphicx}%%
\\begin{document}
    \\input{%s_small}%%
\\end{document}
%}

function [str] = xelatex
str = verbatim;
%{
\\usepackage[no-math]{fontspec}
\\usepackage{xunicode,xltxtra}

%}

function [str] = media9(pdforxelatex)
if strcmp(pdforxelatex, 'xelatex')
    str = '\\usepackage[dvipdfmx]{media9}\n';
else
    str = '\\usepackage[dvipdfmx]{media9}\n';
end

function [str] = movie15
str = verbatim;
%{
\\usepackage{hyperref}%%
\\usepackage[3D]{movie15}%%

%}
