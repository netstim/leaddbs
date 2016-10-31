% writeOFF.m
%
%       $Id$	
%      usage: writeOFF(surf, outName)
%         by: eli merriam
%       date: 10/25/07
%    purpose: 
%
function retval = writeOFF(surf, outName)

if nargin ~= 2
  help writeOFF;
  return
end
% write the OFF format file

% undo what loadSurfOFF did in loading surface into matlab
surf.tris = surf.tris - 1;
surf.tris = cat(1, repmat(3,1,length(surf.tris)), surf.tris', repmat(0,1,length(surf.tris)));

% if the surface has had the vtcs converted from
% world to array coordinates by xfromSurfaceWorld2Array
% then revert to original vertices for saving
if isfield(surf,'originalVtcs')
  surf.vtcs = surf.originalVtcs;
end

fid = fopen(setext(outName,'off'), 'w', 'ieee-be');

fprintf(fid, 'OFF BINARY\n');
fwrite(fid, [surf.Nvtcs surf.Ntris 0], 'int32'); 

% Vertices
fwrite(fid, surf.vtcs', 'float32');

% Faces
fwrite(fid, surf.tris, 'int32');

% Close file
fclose(fid);

