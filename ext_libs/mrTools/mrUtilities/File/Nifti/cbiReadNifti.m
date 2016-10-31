function [data,hdr]=cbiReadNifti(fname,subset,prec,short_nan,verbose)
% [data,hdr]=cbiReadNifti(fname,subset,prec,short_nan,verbose)
% 
% Loads all or part of a Nifti-1 file.
%  Default is to load entire file as a double precision array
%  with dimensions [xSize ySize zSize tSize] where x,y,z,t are
%  defined according to the Nifti standard (x fastest changing 
%  dimension, y second and so on - usually x corresponds to left-right,
%  y to back-front, z to down-up, and t to early-late. 
%  NB: will also work for quasi-5-dimensional data, e.g. 3D arrays of vectors
%  defined on a single time point (t then being the dimensions of the vectors 
%  at each point), i.e. hdr.dim(4)==1 and hdr.dim(5)>1. (0-offset index)
%  Will not work with higher dimensional data sets (hdr.dim(4)>1 & hdr.dim(5)>1).
%  Ignores hdr.dim(6) and higher.
% 
% hdr contains the Nifti header returned by readNiftiHeader.
% 
% Options:
%  subset:    4x1 cell array describing image subset to retrieve. 1-offset (Matlab-style).
%             Examples: 
%             - to retrieve a single z-slice (e.g. 4):  subset={[],[],4,[]}
%             - to retrieve a single voxel time series (e.g. [5,6,7]):  subset={5,6,7,[]}
%             - to retrieve a single volume from a time series (e.g. second volume):  subset={[],[],[],2}
%             - to retrieve a block of voxels from a volume: eg. subset={[4 50],[6 20],[1 10],[]}
%             If data is quasi-5D, subset{4} defines the subset in dim(5)
% 
%  prec:      Storage type of returned data. Legal values are:
%             'native' - data are returned as stored in file (no scaling). Only works for loading
%              contiguous data (entire volumes).
%             'double' - data are returned in Matlab standard double precision format (default)
%             'single' - data are returned in Matlab single precision format
%  short_nan: NaN handling for signed short (int16) data. If 1, will treat -32768 (smallest 
%             representable number) as NaN, reserving -32767..32767 for scaled data; otherwise
%             treats -32786 as a real value. Default is 1. 
% 
% JL 20050223
% $Id$


  
  
if ieNotDefined('subset')
  subset={[],[],[],[]};
end
if ieNotDefined('prec')
  prec=mrGetPref('defaultPrecision');
end
if ieNotDefined('short_nan')
  short_nan=1;
end
if ieNotDefined('verbose')
  verbose=1;
end

% first check for valid filename
if isempty(getext(fname))
  fname = setext(fname,'hdr');
end
if ~isfile(fname)
  disp(sprintf('(cbiReadNifti) UHOH! Could not find file %s',fname));
  data =[];hdr=[];
  return
end

hdr=cbiReadNiftiHeader(fname);

% this is a bit of a hack, but the pixdim
% can be at most 2^15-1 which is a problem
% for base surface images which just
% contains info about each vertex. So
% we check for that condition and fix
% the dimension here
needsRewrite = 0;
maxWidth = 2^15-1;
if (hdr.dim(2) ==1) && (hdr.dim(3) == maxWidth)
  dataDir = dir(sprintf('%s.img',stripext(fname)));
  dataBytes = dataDir(1).bytes/(hdr.bitpix/8);
  if dataBytes ~= hdr.dim(3)
    if verbose, disp(sprintf('(cbiReadNifti) Image dimension 2 greater than nifti max (%i), fixing and rereading %i bytes',hdr.dim(3),dataBytes));end
    hdr.dim(3) = dataBytes;
  end
  % this part is payback for this cheap hack, I did above. Rewrite
  % the file now as a valid nifti image with width = 2^15-1 and
  % enough rows to fit the image, by calling cbiWriteNifti -j.
  needsRewrite = 1;
end

headerdim=hdr.dim(2:5); % Matlab 1-offset - hdr.dim(1) is actually hdr.dim(0)
headerdim(headerdim==0)=1; % Force null dimensions to be 1
is5D=0;
if (hdr.dim(6)>1)
  if (hdr.dim(5)>1)
    error('No support for 5D data with multiple time points!');
  end
  is5D=1;
  headerdim(4)=hdr.dim(6);
  if verbose, disp('5D data set detected');end;
end

loadSize=zeros(4,1);
for n=1:4
  if (length(subset{n})==1)
    subset{n}=[subset{n} subset{n}];    
  elseif (length(subset{n})>2)
    error('subset should be a scalar or 2-vector');
  end
  if (isempty(subset{n}))
    loadSize(n)=headerdim(n);
    subset{n}=[1 loadSize(n)];
  else
    loadSize(n)=subset{n}(2)-subset{n}(1)+1;  
  end
end
if (any(loadSize>headerdim(1:4)))
  error('subset index larger than image dimensions!');
end

% Set read readformat: native or double.
if (strcmp(prec,'native'))
  if (loadSize(1:3)==headerdim(1:3))
    readformat=[hdr.matlab_datatype '=>' hdr.matlab_datatype]; % Loads data in native format
  else
    mrWarnDlg('(cbiReadNifti) native readformat only supported for loading contiguous volumes. Forcing double precision...');
    prec='double';
    readformat=hdr.matlab_datatype;    
  end
else
  readformat=hdr.matlab_datatype;
end

% If native format is unsupported by Matlab, force conversion
switch (hdr.matlab_datatype)
 case 'binary'
  error('Binary format not supported!');
 case 'complex256','float128'
  error('No support for 128-bit data on this platform!');
 case 'complex64'
  if (strcmp(prec,'native'))
    mrWarnDlg('Warning: native data type not supported by Matlab. Forcing conversion to double.');
  end
  readformat='float32';
  prec='double';
 case 'complex128'
  if (strcmp(prec,'native'))
    mrWarnDlg('(cbiReadNifti) Native data type not supported by Matlab. Forcing conversion to double.');
  end
  readformat='float64';
  prec='double';
 case 'RGB'
  error('No support for RGB data.');
 otherwise
  % nothing
end

% Open file and move to beginning of data
fPtr=fopen(hdr.img_name,'r',hdr.endian);
if fPtr==-1
  disp(sprintf('(cbiReadNifti) Could not find file %s',hdr.img_name));
  data=[];hdr=[];return;
end
% If single-file format, offset beyond header
if (hdr.single_file)
  fseek(fPtr,hdr.vox_offset,'bof');
end

dataSize=prod(loadSize);
bytesPerElement=cbiSizeofNifti(readformat);

if (loadSize(1:3)==headerdim(1:3))
  % Loading of contiguous data
  if (loadSize(4)~=headerdim(4))
    % Time series segment - fseek to read origin
    readOrigin=sub2ind(headerdim(1:4)',subset{1}(1),subset{2}(1),subset{3}(1),subset{4}(1))-1;
    fseek(fPtr,readOrigin*bytesPerElement,'cof');
  end
  switch (hdr.matlab_datatype)
   case {'complex64','complex128'}
    [cdata,count]=fread(fPtr,dataSize*2,readformat);
    if (count~=dataSize*2) 
      fclose(fPtr);
      error(['Error reading file ' hdr.img_name]);
    end
    data=complex(cdata(1:2:dataSize*2-1),cdata(2:2:dataSize*2));
    clear cdata;
   otherwise
    % Load entire file
    [data,count]=fread(fPtr,dataSize,readformat);
    if (count~=dataSize) 
      fclose(fPtr);
      mrWarnDlg(sprintf('(cbiReadNifti) File contains %i bytes, but expected %i. Could not read (file corrupted?)',count,dataSize));
      data = [];
      return
    end
  end
  data=reshape(data,loadSize');
else
  % Use fseek to load segments into array 
  volSize=prod(headerdim(1:3));
  % Initialize data array
  if ( prod(loadSize(1:3))>1 & (loadSize(1)<headerdim(1) | loadSize(2)<headerdim(2)))
    data=zeros(loadSize');
    % 1) Non-contiguous data (restrictions in x or y): read a volume at a time and restrict
    readSize=volSize;
    if (strfind(hdr.matlab_datatype,'complex'))
      readSize=readSize*2;
    end
    for t=subset{4}(1):subset{4}(2)
      % load single volume
      [d,count]=fread(fPtr,readSize,readformat);
      if (count~=readSize) 
	fclose(fPtr);
	error(sprintf('(cbiReadNifti) Expected to read %i values from %s, but only read %i',readSize,hdr.img_name,count));
      end
      d=reshape(d,headerdim(1:3)');
      xsub=subset{1}(1):subset{1}(2);
      ysub=subset{2}(1):subset{2}(2);
      zsub=subset{3}(1):subset{3}(2);
      if (strfind(hdr.matlab_datatype,'complex'))
	% Convert real and imaginary parts into complex representation
	cd=complex(d(1:2:readSize/2-1),d(2:2:readSize/2));	
	data(xsub-subset{1}(1)+1,ysub-subset{2}(1)+1,zsub-subset{3}(1)+1,t)=cd(xsub,ysub,zsub);
      else      
	data(xsub-subset{1}(1)+1,ysub-subset{2}(1)+1,zsub-subset{3}(1)+1,t)=d(xsub,ysub,zsub);
      end
    end    
  else
    data=zeros(prod(loadSize),1);
    % 2) Single voxel over time: fseek to voxel index and start time and read a voxel every volOffset
    % 3) Single or multiple slices over time: fseek to slice index and start time and read a block every volOffset
    % Determine file read origin from first voxel index
    readOrigin=sub2ind(headerdim(1:4)',subset{1}(1),subset{2}(1),subset{3}(1),subset{4}(1))-1; % now we're in C-land, hence 0-offset
    % Elements to read every time point
    readSize=prod(loadSize(1:3));
    % Difference between volSize and readSize => offset to seek after every read
    readOffset=volSize-readSize;
    % Current position in data array
    currPos=1;
    if (strfind(hdr.matlab_datatype,'complex'))
      % Every voxel corresponds to two elements in the file
      readOrigin=readOrigin*2;
      readSize=readSize*2;
      readOffset=readOffset*2;
    end
    % Position file at first voxel
    fseek(fPtr,readOrigin*bytesPerElement,'cof');
    for t=subset{4}(1):subset{4}(2)
      % Read one temporal chunk of data
      [d,count]=fread(fPtr,readSize,readformat);
      if (count~=readSize) 
	fclose(fPtr);
	error(['Error reading file ' hdr.img_name]);
      end
      if (readOffset>0)
	% fseek to next time point
	fseek(fPtr,readOffset*bytesPerElement,'cof');
      end
      if (strfind(hdr.matlab_datatype,'complex'))
	% Convert real and imaginary parts into complex representation
	cd=complex(d(1:2:readSize/2-1),d(2:2:readSize/2));
	data(currPos:currPos+readSize/2-1)=cd;
      else
	data(currPos:currPos+readSize-1)=d;
      end
      currPos=currPos+readSize;
    end
    data=reshape(data,loadSize');
  end
end
  
fclose(fPtr);

% Scaling
if (strcmp(prec,'double')||strcmp(prec,'single'))
  if hdr.scl_slope~=0
    % NaN handling for int16 data
    if (hdr.datatype==4 & short_nan) 
      data(data==-32768)=NaN;
    end
    data=data.*hdr.scl_slope+hdr.scl_inter;
  end
  % convert to single precision down here--probably better to read
  % in the data as a single, but that requires changing a lot of
  % code above -jg.
  if strcmp(prec,'single')
    data=single(data);
  end
end

% this part is payback for a cheap hack, I did above. Rewrite
% the file now as a valid nifti image with width = 2^15-1 and
% enough rows to fit the image, by calling cbiWriteNifti
if needsRewrite
  [tf permissions] = isfile(fname);
  if tf && permissions.UserWrite
    if strcmp(questdlg(sprintf('(cbiReadNifti) %s is an old style surface curvature file. Answer yes to update the format of the surface (recommended).',getLastDir(fname)),'Yes','No'),'Yes')
      if verbose, disp(sprintf('(cbiReadNifti) Fixing 1D image dimensions by rewriting %s',fname));end;
      cbiWriteNifti(fname,data,hdr); 
    end
  else
    mrWarnDlg('(cbiReadNifti) This is probably an old style surface file, you are advised to (re)save it to update the format');
  end
end

% 1D curvature file fix: if the width is at maximum, we should reshape it to a 1xn image
if (hdr.dim(3) == maxWidth)
  % check to see if the description has the real width
  [realWidthStart realWidthEnd] = regexp(hdr.descrip,'^[0-9]*:');
  if ~isempty(realWidthStart)
    % get the real width
    realWidth = str2num(hdr.descrip(realWidthStart:realWidthEnd-1));
    % reset the description
    hdr.descrip = hdr.descrip(realWidthEnd+1:end);
    hdr.descrip(end+1:80) = 0;
    % tell user what we are doing
    if verbose, disp(sprintf('(cbiReadNifti) Fixing dimensions of image that was written as %ix%i but actually is 1x%i',hdr.dim(2),hdr.dim(3),realWidth));end;
    % fix the header
    hdr.dim(2) = 1;
    hdr.dim(3) = realWidth;
    % fix the data
    data = data(1:realWidth);
  elseif (hdr.dim(2)==1)
    mrErrorDlg(sprintf('(cbiReadNifti) %s is an old style surface curvature file. You need to re-import the surface',getLastDir(fname)));
  end
end


  
