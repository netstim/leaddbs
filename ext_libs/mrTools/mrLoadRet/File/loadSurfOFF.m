% loadSurfOFF.m
%
%      usage: surf = loadSurfOFF(surffile,<loadOnlyHeader>)
%         by: jonas larsson (adapted to mrLoadRet by eli merriam)
%       date: 09/25/07
%    purpose: Loading an OFF binary surface into Matlab
%   $Id$	
%
% This code was ripped directly from Jonas Larsson's tfiReadOFF code
% the main difference is that this code returns a neat structure,
% rather than a bunch of variables
%
% ***IMPORTANT NOTE*** Jonas' SurfRelax surfaces keep coordinates in "world" coordinates
% that is, they are offset by the offset in the nifti header of the volume
% anatomy that they were created from. To get the indexes into the volume
% anatomy (which is what you usually want), you have to run the function:
%
% surf = xformSurfaceWorld2Array(surf,hdr)
% 
% where hdr is the nifti header of the volume anatomy the surface was generated from
%
% INPUT ARGUMENTS:
% surffile - file name of surface in OFF binary format.
%
% OUTPUT ARGUMENTS:
% vtcs: Nvtcs x 3 matrix of vertex coordinates
% tris: Ntris x 3 matrix of triangle vertex indexes
%       NB!!! Vertex indexes are 1-offset for Matlab compatibility.
%       In the OFF surface vertex indexes are 0-offset.
%
% The surf structure consists of the following variables
% vtcs,tris,Nvtcs,Ntris,Nedges (nParent,nPatch,patch2parent)
%
function surf = loadSurfOFF(surffile,loadOnlyHeader)

% check arguments
if ~any(nargin == [1 2])
  help loadSurfOFF
  return
end

if ieNotDefined('loadOnlyHeader'),loadOnlyHeader = 0;end

% put on .off
surffile = sprintf('%s.off',stripext(surffile));

if ~isfile(surffile)
  surf = [];
  disp(sprintf('(loadSurfOFF) Could not find surface file %s',surffile));
  return
end

surf.filename = surffile;
fid = fopen(surffile, 'r', 'ieee-be');
fl = fgetl(fid);
if (~strcmp(fl,'OFF BINARY'))
  if (strcmp(fl,'#PATCH'))

    % get the surface parent name
    fl = fgetl(fid);
    surf.parentSurfaceName = sscanf(fl,'#parent_surface=%s\n');
    
    % parent surface dimensions
    fl = fgetl(fid); 
    surf.nParent = sscanf(fl,'#parent_dimensions=%i %i %i\n');

    % patch surface dimensions
    fl = fgetl(fid); 
    surf.nPatch = sscanf(fl,'#patch_dimensions=%i %i %i\n');

    % parent vertex indexes tag
    fl = fgetl(fid); 
    
    % read patch info
    [a,c] = fscanf(fid,'#%i %i\n',[2 surf.nPatch(1)]);
    if (c ~= surf.nPatch(1)*2)
      disp(c)
      mrErrorDlg(sprintf('(loadSurfOFF) Error reading file %s',surffile));
    end
    
    % matlab 1-offset
    surf.patch2parent = a'+1; 
    while (~strcmp(fl,'OFF BINARY'))
      fl= fgetl(fid);
    end
  else
    disp('WARNING!! Magic number (OFF BINARY) not found - this file may not be in the correct format!');
  end
end

[Ninfo, count] = fread(fid, 3, 'int32');
if (count~=3), mrErrorDlg(sprintf('(loadSurfOFF) Error reading file %s',surffile)); end
surf.Nvtcs  = Ninfo(1);
surf.Ntris  = Ninfo(2);
surf.Nedges = Ninfo(3);

% return here if we only want the header
if loadOnlyHeader
  fclose(fid);
  return
end

[surf.vtcs, count] = fread(fid, surf.Nvtcs*3, 'float32');
if (count ~= 3*surf.Nvtcs), mrErrorDlg(sprintf('(loadSurfOFF) Error reading file %s',surffile)); end

[surf.tris, count] =  fread(fid, surf.Ntris*5, 'int32');
if (count ~= 5*surf.Ntris), mrErrorDlg(sprintf('(loadSurfOFF) Error reading file %s',surffile)); end

surf.vtcs = reshape(surf.vtcs, 3, surf.Nvtcs);
surf.tris = reshape(surf.tris, 5, surf.Ntris);

surf.vtcs = surf.vtcs';

% first entry is # of vertices per face (always 3); last entry is # of colors per face (only 0 allowed)
surf.tris = surf.tris(2:4,:); 
% 1-offset for matlab
surf.tris = surf.tris'+1; 

fclose(fid);

return;

