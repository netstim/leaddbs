function hdr=cbiSetNiftiQform( hdr, mtx44, trns )
% hdr=cbiSetNiftiQform( hdr, mtx44 )
%  where mtx44 is a 4x4 homogeneous matrix
% - OR -
% hdr=cbiNiftiSetQform( hdr, R, T )
%  where R is a 3x3 rotation matrix and T a translation 3-vector
% 
% Sets qform44 and updates quaternions, pixdim, and qoffsets accordingly.
% Sets qcode=1

if (nargin<2)
  error('Must specify a header and a qform matrix (or a rotation and translation matrix')
end
  
if (nargin==3)
  m=eye(4);
  if (size(mtx44)~=[3 3])
    error ('rotation matrix must be 3x3');
  end
  m(1:3,1:3)=mtx44;
  if (length(trns)~=3)
    error('translation vector must be a 3-vector');
  end
  m(1:3,4)=trns(:);
  mtx44=m;
end

hdr.qform44=mtx44;
[quatern, qfac, pixdim]=cbiHomogeneousToQuaternion(hdr.qform44);

hdr.quatern_b=quatern(1);
hdr.quatern_c=quatern(2);
hdr.quatern_d=quatern(3);

hdr.pixdim(1)=qfac;
hdr.pixdim(2:4)=pixdim;

hdr.qoffset_x=hdr.qform44(1,4);
hdr.qoffset_y=hdr.qform44(2,4);
hdr.qoffset_z=hdr.qform44(3,4);

hdr.qform_code=1;
