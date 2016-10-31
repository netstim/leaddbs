function varargout=cbiSwapNiftiDimensions(varargin);
% [newdata,newhdr,swapMatrix]=cbiSwapNiftiDimensions(data,hdr,swapvect);
% - OR -
% swapMatrix=cbiSwapNiftiDimensions(inputfile,outputfile,swapvect,[outputdatatype]);
% 
% Swaps the dimensions of a Nifti file or data set.
% Modifies the qform44 and quaternions accordingly.
% If sform44 is non-identity, also modifies the sform.
% Equivalent to avwswapdim except it (hopefully correctly) modifies qform44/sform44 matrices
% NOTE: Only modifies the first 3 dimensions of a file - but works with data with 4 or more dimensions.
%
% swapvect: a 3-vector with entries +/-(1,2,3) indicating axis changes/flips, e.g.:
% - flip x (dimension 1):      swapvect=[-1 2 3]
% - exchange x and z:          swapvect=[3 2 1]
% - exchange x with flipped y: swapvect=[-2 1 3]
%
% NOTE: if the image is acquired in one of the cardinal directions (axial, sagittal, or coronal 
% slice prescription) swapvect can be omitted, in which case the default is to swap the image axes 
% to conform to a 'standard' right-handed coordinate system as follows:
% - dimension i (fastest changing dimension) corresponds to x (increasing from left to right)
% - dimension j (second fastest changing dimension) corresponds to y (increasing from back to front)
% - dimension k (third fastest changing dimension) corresponds to z (increasing from bottom to top)
% 
% Optionally returns swap matrix
%
  
% File or data input?

isfile=0;
isarray=0;
if (isstr(varargin{1}) & isstr(varargin{2}))
  % Input is a file
  isfile=1;
  [data,hdr]=cbiReadNifti(varargin{1});
  datasize=size(data);
else
  % Input is a Matlab data array
  if (nargout<2) 
    error('Not enough outputs')
  end
  if (~isstruct(varargin{2}))
    error('Second input must be a header struct')
  end
  hdr=varargin{2};
  isarray=1;
  datasize=size(varargin{1});
end

if (nargin>=3)
  swapvect=varargin{3};
else
  % check that image is in cardinal planes
  image_axes=hdr.qform44(1:3,1:3);
  ndims_per_axis=sum(abs(image_axes)>0);
  if (any(ndims_per_axis>1))
    error('Image planes not in cardinal planes - image must be rotated to canonical orientation. Aborting...')
  else
    swapvect=[0 0 0];
    swapmat=inv(image_axes);
    for n=1:3
      % row index of nonzero entry in each column == axis index with sign
      ri=find(swapmat(:,n));
      swapvect(n)=ri*sign(swapmat(ri,n));
    end
  end
  disp(['Extracted permutation vector (swapvect): ' num2str(swapvect)])
end

if (length(swapvect)~=3)
  error('Permutation (swap) vector must have 3 dimensions!');
end

if (~find(abs(swapvect)==1) | ~find(abs(swapvect)==2) | ~find(abs(swapvect)==3 ))
  error('swapvect MUST contain 1, 2, and 3')
end

% Generate swap matrix P
% Old coordinates: X
% New coordinates after flipping: Y = P*X
% Scanner coordinates: S = Q*X where Q is the qform44
% Hence the new qform R is given by
% S = Q*X = Q*inv(P)*Y 
% R = Q*inv(P)
% and similarly, the new sform Z is given by
% Z = W*inv(P) where W is the old sform
% NOTE: swapMatrix == inv(P)

swapMatrix=eye(4);
id44=eye(4);
for n=1:3
  % axis changes correspond to exchanging columns
  srccol=swapvect(n);
  sgn=sign(srccol);
  srccol=abs(srccol);
  swapMatrix(:,n)=sgn*id44(:,srccol);
  if (sgn<0)
    % axis inversions correspond to inverting sign and adding axis dimension-1 
    % (because this matrix operates on array coordinates which are zero-offset integers)
    % (i.e. voxel i==2 when flipped becomes voxel i==datasize-1-2)
    swapMatrix(srccol,4)=datasize(srccol)-1;
  end
end
disp(swapMatrix)

% Calculate new qform & sform. Qform is needed for pixdim etc.
if (~isfield(hdr,'qform44') | isempty(hdr.qform44))
  hdr.qform44=eye(4);
end

newqform44=hdr.qform44*swapMatrix;
% Set new qform and update quaternions, qoffset, and pixdim
fliphdr=cbiSetNiftiQform( hdr, newqform44 );

if (isfield(hdr,'sform44') & ~isequal(eye(4),hdr.sform44))
  newsform44=hdr.sform44*swapMatrix;
  fliphdr=cbiSetNiftiSform( fliphdr, newsform44 );
end

% Flip data
permutevect=1:length(datasize);
permutevect(1:3)=abs(swapvect);
if (isarray)
  flipdata=permute(varargin{1},permutevect);
else
  flipdata=permute(data,permutevect);  
end
for n=1:3
  if (swapvect(n)<0)
    flipdata=flipdim(flipdata,n);
  end
end

% Update flipped header
fliphdr=cbiCreateNiftiHeader(fliphdr,flipdata);

if (isarray)
  varargout{1}=flipdata;
  varargout{2}=fliphdr;
  if (nargout==3) 
    varargout{3}=swapMatrix;
  end
else
  % Save if desired
  if (nargin==4)
    [b,h]=cbiWriteNifti(varargin{2},flipdata,fliphdr,varargin{4});    
  else
    [b,h]=cbiWriteNifti(varargin{2},flipdata,fliphdr);
  end
  if (nargout==1) 
    varargout{1}=swapMatrix;
  end
end


