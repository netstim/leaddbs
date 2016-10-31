function xform=mlrXFormFromHeader( filename, convtype )
% xform=mlrXFormFromHeader( filename, convtype )
%
% Returns a transformation matrix relating array (Matlab matrix)
% coordinates to Nifti or SurfRelax 'world' coordinates.
% filename : Analyze or Nifti image or header file
% convtype : conversion to extract, one of:
%            'array2world', 'a2w' : Matlab array coordinates to world (SurfRelax) coordinates
%            'array2nifti', 'a2n', 'array2qform', 'a2q' : Matlab array coordinates to Nifti qform coordinates
%            'world2array', 'w2a' : World (SurfRelax) coordinates to Matlab array coordinates
%            'nifti2array', 'n2a', 'qform2array', 'q2a' : Nifti qform coordinates to Matlab array coordinates
%            'array2sform', 'a2s' : Matlab array coordinates to Nifti sform coordinates
%            'sform2array', 's2a' : Matlab array coordinates to Nifti sform coordinates
%
% In all cases, the returned xform maps a homogeneous coordinate in
% column vector format Pfrom to Pto as follows:
%
% Pto = xform * Pfrom
% 
% i.e. to get the Nifti coordinate of Pfrom=[1 2 4], compute
% 
% Pto = xform * ([Pfrom 1]')
%
% Note: Surfaces and volumes made with SurfRelax version<2 do not use
% Nifti coordinates. This will change in next release (v.2), making the
% world coordinates obsolete.
%
% Note, that one change was made from Jonas' original code, in that
% filename can be a passed in hdr -jlg
% 

if (nargin<2)
  help(mfilename)
  return
end

if isstruct(filename)
  hdr = filename;
else
  hdr=cbiReadNiftiHeader(filename);
end


switch (convtype)
 case {'array2world', 'a2w'}
  % calls local function below
  xform = array2world(hdr);
  
 case {'array2nifti', 'a2n', 'array2qform', 'a2q'}
  xform = hdr.qform44;
  % Add -1 for Matlab 1-offset
  xform(1:3,4)=xform(1:3,4) - 1;
 
 case {'world2array', 'w2a'}
  xform = inv(array2world(hdr));
  
 case {'nifti2array', 'n2a', 'qform2array', 'q2a'}
  xform = hdr.qform44;
  % Add -1 for Matlab 1-offset
  xform(1:3,4)=xform(1:3,4) - 1;
  xform=inv(xform);
  
 case {'array2sform', 'a2s'}
  xform = hdr.sform44;
  % Add -1 for Matlab 1-offset
  xform(1:3,4)=xform(1:3,4) - 1;

 case {'sform2array', 's2a'}
  xform = hdr.sform44;
  % Add -1 for Matlab 1-offset
  xform(1:3,4)=xform(1:3,4) - 1;
  xform = inv(xform);
  
 otherwise
  help(mfilename)
  error('Unknown conversion.')
end

return

function xform=array2world(hdr)
% Old Surfrelax world coordinates are computed as  x_array*xAspect_ - xOrigin_ 
% or as (x_world+xOrigin_)/xAspect_); }; 
% where aspect is pixel dimension and origin is pixdim[5-7]
% (C-style, ie pixdim(6:8) matlab)

xform=eye(4);
xform(1:3,1:3)=diag(hdr.pixdim(2:4));
% Note -1 for Matlab 1-offset
xform(1:3,4)=-hdr.pixdim(6:8) - 1;


return
