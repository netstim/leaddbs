% getfid.m
%
%      usage: getfid(fidname,<verbose=1>,<zeropad=0>,<movepro=0>,<kspace=0>,<swapReceiversAndSlices=1>)
%         by: justin gardner
%       date: 05/08/03
%    purpose: reads k-space data from fid file and transforms into an image
%             and transforms it into image. Dimensions are
%             set to x,y,slice,vol,coil (note that vol and coil
%             dimensions are swapped relative to what getfidk returns)
% 
%             optional arguments
%             verbose=0 to 1 for verbose info
%             zeropad=0 set to the size you want to zeropad out to (e.g. 256)
%             movepro=0 set to how much you want to movepro (default=0)
%             kspace=0, set to 1 if you want k-space rather than image data
%             swapReceiversAndSlices=1, swaps receivers and slices to keep%
%               consistent dimension ordering. For 2D images with multiple receivers, Varian
%               encodes the image with the slices and receivers swapped. This
%                program (but not getfidk) fixes that swap with this argument set.
%
function d = getfid(fidname,varargin)

% check input arguments
if nargin<1
  help getfid;
  return
end

verbose=[];zeropad=[];movepro=[];kspace=[];movepss=[];
oldArgNames = {'verbose','zeropad','movepro','kspace'};
% check for numeric arguments (old way of calling
for i = 1:length(varargin)
  if isnumeric(varargin{i})
    % if it is numeric then change it into something that getArgs can parse
    if isempty(varargin{i})
      varargin{i} = sprintf('%s=[]',oldArgNames{i});
    else
      varargin{i} = sprintf('%s=%s',oldArgNames{i},mynum2str(varargin{i}));
    end
  else 
    % if non-numeric then stop doing this, since we are using new style getArgs
    break
  end
end
% now run getArgs
getArgs(varargin,{'verbose=0','zeropad=0','movepro=0','kspace=0','swapReceiversAndSlices=1','movepss=0'});

% read the k-space data from the fid
if (verbose),disppercent(-inf,sprintf('(getfid) Reading %s...',fidname));end
d = getfidk(fidname,'verbose',verbose);
if (verbose),disppercent(inf,sprintf('done.\n',fidname));end
% if it is empty then something has failed
if (isempty(d.data))
  return
end

d.dim = size(d.data);

% get fidinfo
[xform info] = fid2xform(fidname);

if info.receiversAndSlicesSwapped && swapReceiversAndSlices
  if size(d.data,4) > 1
    if verbose, disp(sprintf('(getfid) Swapping receivers and slices'));end
    d.data = reshape(d.data,size(d.data,1),size(d.data,2),size(d.data,4),size(d.data,3),size(d.data,5));
    d.data = permute(d.data,[1 2 4 3 5]);
  end
end

% find the needed phase shift to add to the image if we need to movepro
if movepro ~= 0
  % figure out how much we have to shift
  proshift = movepro/(info.procpar.lro/d.dim(1));
  % and create the shift in 2D
  phaseshift2 = proshift*(0:2*pi./d.dim(1):2*pi)';
  phaseshift2 = phaseshift2(2:end);
  phaseshift2 = phaseshift2*ones(1,d.dim(2));
  phaseshift2 =  exp(sqrt(-1)*phaseshift2);
  % and create the shift in 3D
  phaseshift3 = proshift*(0:2*pi./d.dim(1):2*pi)';
  phaseshift3 = phaseshift3(2:end);
  phaseshift3 = repmat(phaseshift3,[1 d.dim(2) d.dim(3)]);
  phaseshift3 = exp(sqrt(-1)*phaseshift3);
  if (verbose),disp(sprintf('(getfid) Shifting pro by: %f',movepro));end
else
  phaseshift2 = 1;
  phaseshift3 = 1;
end

% see if we have to move the pss (only available for 3D volumes)
if movepss ~= 0
  if ~info.acq3d || info.fftw3dexe_processed
    disp(sprintf('(getfid) !!! Unable to move pss for 2D data. Ignoring desired pss shift of: %f !!!',movepss));
  else
    % and create the shift in 3D
    pssshift = movepss/(info.procpar.lpe2/d.dim(3));
    pssPhaseshift3 = pssshift*(0:2*pi./d.dim(3):2*pi)';
    pssPhaseshift3 = pssPhaseshift3(2:end);
    pssPhaseshift3 = repmat(pssPhaseshift3,[1 d.dim(1) d.dim(2)]);
    pssPhaseshift3 = permute(pssPhaseshift3,[2 3 1]);
    pssPhaseshift3 = exp(sqrt(-1)*pssPhaseshift3);
    % composite with phaseshift from movepro
    phaseshift3 = phaseshift3.*pssPhaseshift3;
  end
end

% everything is ok, then transform data
if(verbose),disppercent(-inf,'(getfid) Transforming data');end

% preallocate space for data
if zeropad
  data = nan(zeropad,zeropad,size(d.data,3),size(d.data,5),size(d.data,4));
else
  data = nan(size(d.data,1),size(d.data,2),size(d.data,3),size(d.data,5),size(d.data,4));
end

% decide whether we need to do a 3D transform or 2D
if info.acq3d && ~info.fftw3dexe_processed
  % 3D FFT
  for j = 1:size(d.data,4)
    for k = 1:size(d.data,5)
      data(:,:,:,k,j) = myfft3(d.data(:,:,:,j,k).*phaseshift3,kspace);
    end
  end
else
  % 2D FFT
  for i = 1:size(d.data,3)
    for j = 1:size(d.data,4)
      for k = 1:size(d.data,5)
	% need to zeropad or not
	if zeropad
	  % zeropad the data
	  thisdata=zeros(zeropad,zeropad);
	  % get this image (apply phaseshift for shifting pro - if set to 0 then
	  % phaseshift is just 1 and doesn't do anything)
	  thisdata(1:d.dim(1),1:d.dim(2)) = squeeze(d.data(:,:,i,j,k)).*phaseshift2;
	  % fft 
	  data(:,:,i,k,j) = myfft(thisdata,kspace);
	else
	  % get this image (apply phaseshift for shifting pro - if set to 0 then
	  % phaseshift is just 1 and doesn't do anything)
	  data(:,:,i,k,j) = myfft(d.data(:,:,i,j,k).*phaseshift2,kspace);
	end
      end
      % percent done
      if (verbose) disppercent(calcPercentDone(i,size(d.data,3),j,size(d.data,4)));end
    end
  end
end

d.data = data;
d.dim = size(data);

% if zeropad fix up some stuff
if zeropad
  d.dim(1) = zeropad;d.dim(2) = zeropad;
end
d.zeropad = zeropad;

d.info = info;

if (verbose), disppercent(inf); end

%%%%%%%%%%%%%%%
%    myfft    %
%%%%%%%%%%%%%%%
function data = myfft(data,kspace)

% this just simply takes the 2D fft, shifts and gets the real part of the data.
% if kspace is set, it does nothing
if ~kspace
  data = fftshift(abs(fft2(data)))/(size(data,1)*size(data,2));
end

%%%%%%%%%%%%%%%%
%    myfft3    %
%%%%%%%%%%%%%%%%
function data = myfft3(data,kspace)

% this just simply takes the 3D fft, shifts and gets the real part of the data.
% if kspace is set, it does nothing
if ~kspace
  data = fftshift(abs(fftn(data)))/(size(data,1)*size(data,2)*size(data,3));
end