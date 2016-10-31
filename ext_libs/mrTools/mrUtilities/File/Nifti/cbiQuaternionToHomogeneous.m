function varargout=cbiQuaternionToHomogeneous(varargin)
% [M,R,T]=cbiQuaternionToHomogeneous( quatern, pixdim, qfac, qoffset )
% - OR -
% [M,R,T]=cbiQuaternionToHomogeneous( hdr )
%
% Computes the 4x4 homogeneous transformation matrix from quaternion coefficients quatern, 
% voxel dimensions pixdim, coordinate system orientation factor qfac and 
% translations qoffset. The format of these is as follows:
%   quatern: 3x1-vector [b,c,d]; coefficient a is calculated as a = +sqrt(1.0-b*b+c*c+d*d).
%            (i.e. assumes that quaternion is normalized).
%   pixdim:  x,y,z voxel dimensions (corresponding to NIFTI hdr.pixdims(2:4) [1-offset])
%   qfac:    1 or -1 (corresponding to NIFTI hdr.pixdims(1) [1-offset])
%   qoffset: 3x1-vector with x,y,z offsets (corresponding to NIFTI [hdr.qoffset_x;hdr.qoffset_y;hdr.qoffset_z])
% 
% Alternatively, extracts these fields from a NIFTI-compliant 
% header struct hdr. 
% 
% If 3 output arguments are given, also returns the rotation matrix R and translation vector T of M.
% Jonas Larsson 2005-02-28
  
% 
% The quaternion (a,b,c,d) is assumed normalized: a = +sqrt(1.0-b*b+c*c+d*d).
% The (b,c,d) values are stored in the (quatern_b,quatern_c,quatern_d) 
% fields, qfac is stored in pixdims[0]
%
%       [ a*a+b*b-c*c-d*d   2*b*c-2*a*d       2*b*d+2*a*c     ] [1]
%   R = [ 2*b*c+2*a*d       a*a+c*c-b*b-d*d   2*c*d-2*a*b     ] [1]
%       [ 2*b*d-2*a*c       2*c*d+2*a*b       a*a+d*d-c*c-b*b ] [qfac]

  if (nargin==1)
    % extract from header
    hdr=varargin{1};
    if ~isstruct(hdr)
      error('Must provide a valid NIFTI-1 header struct!');
    end
    qfac = hdr.pixdim(1);
    b = hdr.quatern_b;
    c = hdr.quatern_c;
    d = hdr.quatern_d;
    pixdim=hdr.pixdim(2:4); % 1-offset in matlab, so pixdims[0] becomes pixdims(1)
    qoffset=[hdr.qoffset_x;hdr.qoffset_y;hdr.qoffset_z];
  elseif (nargin==4)
    b=varargin{1}(1);
    c=varargin{1}(2);
    d=varargin{1}(3);    
    pixdim=varargin{2};
    qfac=varargin{3};
    qoffset=varargin{4};
  else
    error('Wrong number of input arguments');
  end

  if (abs(qfac)<1e-10)
    % This used to make a warning, but the nifti documentation
    % suggests that it is correct to assume that qfac should
    % be set to 1 if we are reading NIFTI-1 headers...so,
    % the warning has been removed:
    %
    %Method 2 uses a factor qfac which is either -1 or 1; qfac is
    %stored in the otherwise unused pixdim[0].  If pixdim[0]=0.0 (which
    %should not occur), we take qfac=1.  Of course, pixdim[0] is only used
    %when reading a NIFTI-1 header, not when reading an ANALYZE 7.5
    %header.
%    disp(['(cbiQuaternionToHomgenous) Invalid qfac (' num2str(qfac) '). Assuming qfac==1']);
    qfac=1;
  end

  a=sqrt(1.0-(b*b+c*c+d*d));

  if (any(imag(a)))
%    disp('(cbiQuaternionToHomogenous) Imaginary residual found.')
    %disp(['Quaternion component a=' num2str(a)])
    %disp('Assuming this is roundoff error, and forcing a=0.0.');
    %disp('If this is not what you want, you may want to reconstruct the quaternion manually from the header.')
    a=0;
  end
  
  R = [ a*a+b*b-c*c-d*d   2*(b*c-a*d)       2*(b*d+a*c); ...
	2*(b*c+a*d)       a*a+c*c-b*b-d*d   2*(c*d-a*b); ...
	2*(b*d-a*c)       2*(c*d+a*b)       a*a+d*d-c*c-b*b ];

  
  M=eye(4);
  M(1:3,1:3) = (R * diag([pixdim(1) pixdim(2) qfac*pixdim(3)]));
  T = qoffset(:);
  M(1:3,4) = qoffset(:);

  varargout{1}=M;
  
  if (nargout==3)
    varargout{2}=R;
    varargout{3}=T;
  end
