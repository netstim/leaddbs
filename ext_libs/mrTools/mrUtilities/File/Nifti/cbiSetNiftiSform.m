function hdr=cbiSetNiftiSform( hdr, mtx44, trns )
% hdr=cbiSetNiftiSform( hdr, mtx44 )
%  where mtx44 is a 4x4 homogeneous matrix
% - OR -
% hdr=cbiSetNiftiSform( hdr, R, T )
%  where R is a 3x3 rotation matrix and T a translation 3-vector
% 
% Sets sform44 and updates srow_x/y/z. 
% Sets scode=1 

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

hdr.sform44=mtx44;

hdr.srow_x=hdr.sform44(1,:);
hdr.srow_y=hdr.sform44(2,:);
hdr.srow_z=hdr.sform44(3,:);

% only reset sform_code if had been 0 (eg unset)
% otherwise don't change it (because may be =3, not =1)
if isequal(hdr.sform_code,0) || isempty(hdr.sform_code)
  hdr.sform_code=1;
end

