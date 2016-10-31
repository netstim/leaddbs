function [curv, fnum] = freesurfer_read_curv(fname)

% freesurfer_read_curv - FreeSurfer I/O function to read a curvature file
%
% [curv, fnum] = freesurfer_read_curv(fname)
% 
% reads a binary curvature file into a vector
%
% After reading an associated surface, with freesurfer_read_surf, try:
% patch('vertices',vert,'faces',face,...
%       'facecolor','interp','edgecolor','none',...
%       'FaceVertexCData',curv); light
% 
% See also freesurfer_write_curv, freesurfer_read_surf, freesurfer_read_wfile


% $Revision$ $Date$

% Copyright (C) 2000  Darren L. Weber
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.

% History:  08/2000, Developed at MGH, Boston
%           01/2004, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ver = '$Revision$ $Date$';
fprintf('FREESURFER_READ_CURV [v %s]\n',ver(11:15));

% open it as a big-endian file
fid = fopen(fname, 'rb', 'b');
if (fid < 0),
    str = sprintf('could not open file: %s', fname);
    error(str);
end

fprintf('...reading surface file: %s\n', fname);
tic;

vnum = freesurfer_fread3(fid);
NEW_VERSION_MAGIC_NUMBER = 16777215;
if (vnum == NEW_VERSION_MAGIC_NUMBER),
    fprintf('...reading new version (float)\n');
    vnum = fread(fid, 1, 'int32');
    fnum = fread(fid, 1, 'int32');
    vals_per_vertex = fread(fid, 1, 'int32');
    curv = fread(fid, vnum, 'float');
else
    fprintf('...reading old version (int16)\n');
    fnum = freesurfer_fread3(fid);
    curv = fread(fid, vnum, 'int16') ./ 100;
end

fclose(fid);
t=toc; fprintf('...done (%6.2f sec)\n\n',t);

return
