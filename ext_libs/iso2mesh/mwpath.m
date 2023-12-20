function tempname=mwpath(fname)
%
% tempname=meshtemppath(fname)
%
% get full temp-file name by prepend working-directory and current session name
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input:
%    fname: input, a file name string
%
% output:
%    tempname: output, full file name located in the working directory
%
%    if global variable ISO2MESH_TEMP is set in 'base', it will use it
%    as the working directory; otherwise, will use matlab function tempdir
%    to return a working directory.
%
%    if global variable ISO2MESH_SESSION is set in 'base', it will be
%    prepended for each file name, otherwise, use supplied file name.
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

temp = getvarfrom({'caller','base'}, 'ISO2MESH_TEMP');
session = getvarfrom({'caller','base'}, 'ISO2MESH_SESSION');

if isempty(temp)
    if isunix
        username = getenv('USER'); % for Linux/Unix/Mac OS
    else
        username = getenv('USERNAME'); % for windows
    end
	iso2meshdir = fullfile(tempdir, ['iso2mesh-' username]);
else
    iso2meshdir = fullfile(temp, session);
end

if ~isfolder(iso2meshdir)
    mkdir(iso2meshdir);
end

if(nargin==0)
    tempname = iso2meshdir;
else
    tempname = fullfile(iso2meshdir, fname);
end
