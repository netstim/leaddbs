function [] = idtf2u3d(fn)
%IDTF2U3D   Convert IDTF to U3D file.
%
% usage
%   IDTF2U3D
%   IDTF2U3D(IDTF_filename)
%   IDTF2U3D(IDTF_filename, U3D_filename)
%
% optional input
%   idtffile = filename string for IDTF file (default = 'matfig.idtf')
%   u3dfile = filename string for U3D file (default = 'matfig.u3d')
%
% output
%   Conerts the IDTF file into a U3D file which is saved to the disk.
%
% note
%   If only IDTF file name is provided without extension, then the '.idtf'
%   file extension is appended and the U3D file uses the same name with the
%   '.u3d' file extension appended.
%
%   If only the IDTF file name is provided with the extension '.idtf', then
%   the U3D file uses the same name with the '.idtf' extension replaced by
%   the extension '.u3d'.
%
%   If both file names are provided, any one without the appropriate
%   extension gets appended with that extension ('.idtf' and '.u3d',
%   respectively).
%
% See also FIG2U3D, FIG2PDF3D, FIG2IDTF.
%
% Based on MESH_TO_LATEX.m by Alexandre Gramfort, which is part of
% "Matlab mesh to PDF with 3D interactive object"
% which is Copyright (c) by Alexandre Gramfort under the BSD License
% The link on the MATLAB Central File Exchange is:
% http://www.mathworks.com/matlabcentral/fileexchange/25383-matlab-mesh-to-pdf-with-3d-interactive-object

% depends
%   clear_file_extension, check_file_extension
%   IDTFConverter executables in ./bin directory

%% input
if nargin < 1
    idtffile = 'matfig.idtf';
    u3dfile = 'matfig.u3d';
end

%% filenames & extensions
idtffile = [fn, '.idtf'];
u3dfile = [fn, '.u3d'];

%% prepare command
execdir = [fileparts(mfilename('fullpath')), filesep, 'bin',filesep, computer, filesep];
ea_libs_helper(execdir);

if ispc
    IDTF = [execdir, 'IDTFConverter.exe'];
else
    IDTF = [execdir, 'IDTFConverter'];
end

%% idtf -> u3d conversion
cmd = [IDTF, ' -input ', ea_path_helper(idtffile), ' -output ', ea_path_helper(u3dfile)];
disp(cmd);

[status, result] = system(cmd);

disp(result)
if status ~= 0
    warning('idtf2u3d:conversion',...
          'IDTFConverter executable returned with error.')
end
