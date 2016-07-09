function [] = u3d_in_latex(fname, media9_or_movie15)
%U3D_IN_LATEX   LaTeX code which includes a U3D file.
%
% usage
%   U3D_IN_LATEX(fname, media9_or_movie15)
%
% input
%   fname = LaTeX file name (a '_small.tex' will be appended)
%   media9_or_movie15 = LaTeX package to use for including the 3D file
%
% output
%   Saves a LaTeX file including the U3D file
%
% File:      u3d_in_latex.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.21
% Language:  MATLAB R2012a
% Purpose:   produce LaTeX code including a U3D file
%
% acknowledgment
%   Based on write_latex by Alexandre Gramfort.
%   This can be found on the MATLAb Central File Exchange:
%       http://www.mathworks.com/matlabcentral/fileexchange/25383-matlab-mesh-to-pdf-with-3d-interactive-object
%   and is covered by the BSD License.

fname = clear_file_extension(fname, '.tex');
texfile = [fname, '_small.tex'];
disp(['Writing: ', texfile] );

%[~, u3dfname] = fileparts(fname);
if strcmp(media9_or_movie15, 'media9')
    content = latex_include_using_media9;
elseif strcmp(media9_or_movie15, 'movie15')
    content = latex_include_using_movie15;
else
    error('media9 or movie15 are the only available choices.')
end

fid = fopen(texfile, 'w');
    fprintf(fid, content, fname, fname, fname);
fclose(fid);

%fprintf(fid, content, fname, fname, ceil(va),...
%        up(1), up(2), up(3), ro, fname);

function [str] = latex_include_using_media9
% VWS file name (w/o extension) = string,
% 2D substitute image file name (w/o extension) = string,
% U3D file name (w/o extension) = string
str = verbatim;
%{
\\includemedia[%%
width=\\linewidth,%%
height=\\linewidth,%%
activate=pagevisible,%%
deactivate=pageinvisible,%%
3Dviews=%s.vws,%%
3Dtoolbar]{%%
    \\includegraphics[width=\\linewidth]{%s.png}%%
}{%s.u3d}%%
%}

function [str] = latex_include_using_movie15
% u3dfile = string,
% u3dfile = string,
% u3dfile = string
str = verbatim;
%{
\\includemovie[%%
poster=%s.png,%%
3Dviews2=%s.vws,%%
autoplay,%%
toolbar]{\\linewidth}{\\linewidth}{%s.u3d}
%}
