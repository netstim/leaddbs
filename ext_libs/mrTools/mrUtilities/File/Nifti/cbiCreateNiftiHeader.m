function hdr=cbiCreateNiftiHeader(varargin)
% hdr=cbiCreateNiftiHeader(varargin)
%
% Generates a NIFTI-compatible Matlab header struct that can be saved with cbiWriteNiftiHeader.
% 
% SYNTAX:
% hdr=cbiCreateNiftiHeader
%     No input arguments: generates a bare-bone header struct with correctly sized fields. 
%     Sets endianness and default data type (float), with voxel size=1 mm isotropic.     
% hdr=cbiCreateNiftiHeader(data)
%     Input is a data array: generates a header to match data (essentially as option (1) but also sets hdr.dims)
% hdr=cbiCreateNiftiHeader(hdr)
%     Input is a header struct: adds missing fields to an existing header and ensures consistency.
% hdr=cbiCreateNiftiHeader(hdr_field1,value1,hdr_field2,value2,...)
%     Input is a set of parameter-value pairs, e.g. ( 'srow_x', [0 0 0 1] ). Sets the corresponding field.
%
%     The output is checked to make sure it is consistent (within limits). Specifically:
%     - If a data array is given as input, hdr.dim will be set from this (overriding parameter-value pairs)
%     - Quaternions will be set as follows:
%       a) if neither qform44, or quatern_[b/c/d] are specified, will use existing hdr info 
%          (if valid, qform44 takes precedence over quaternions)
%       b) specifying a quaternion/qoffset as a parameter sets qform44 accordingly
%       c) specifying a qform44 sets quaternions accordingly (overriding any quaternion parameters)
%       d) if qform44 is empty, qform_code=0, else qform_code=1
%     - Sform data will be set as follows:
%     - a) if neither sform44 or srow_[x/y/z] are specified, will use  existing hdr info 
%          (if valid, sform44 takes precedence over srow_[x/y/z]
%     - b) specifying a srow_ paramater sets sform44 accordingly
%     - c) specifying sform44 sets srows accordingly (overriding any srow parameters)
%     - d) if sform44 is empty, sform_code=0, else sform_code is left unchanged
%     - Setting the matlab_datatype also sets the NIFTI data type and bitpix appropriately (and vice versa). 
%       (note: bitpix is determined by the datatype; setting it directly has no effect) 
%     - Setting the file name also sets the magic number
% 
% Options 2-4 can be combined. The data and header input arguments must come before the parameter-value pairs.
  
  if (nargin==0)
    % create empty header
    hdr=generateEmptyNiftiHeader;
    return
  end
  
  % Check if the input contains a data or existing header struct
  input_hdr=0;
  input_data=0;
  for currarg=1:min([nargin;2])
    if (isstruct(varargin{currarg}))      
      input_hdr=currarg;
    elseif (~isstr(varargin{currarg}))
      if (currarg==1)
	% must be a data array - parameter/value pairs always start with a string
	input_data=currarg;
      else
	% can only be a data array if the first argument was the header
	if (input_hdr==1)
	  input_data=currarg;
	end
      end
    end
  end

  % Make a header
  if (input_hdr)
    hdr=generateEmptyNiftiHeader(varargin{input_hdr});
  else
    hdr=generateEmptyNiftiHeader;
  end

  use_quatern=0;
  use_qform=0;
  use_srow=0;
  use_sform=0;
  
  % Set parameter-value pairs
  argstart=max([input_hdr input_data])+1;
  if (argstart<nargin)
    if (rem(nargin-argstart+1,2))
      error('Missing value for parameter-value pairs');
    end
    
    for currarg=argstart:2:nargin
      if (~isstr(varargin{currarg}))
	error('Incorrect parameter specification (expected a field name)');
      end
      hdr=setfield(hdr,varargin{currarg},varargin{currarg+1});
      if ~isempty(findstr('srow',varargin{currarg}))
	use_srow=1;
      end
      if ~isempty(findstr('sform44',varargin{currarg}))
	use_sform=1;
      end
      if (~isempty(findstr('quatern',varargin{currarg})) | ~isempty(findstr('qoffset',varargin{currarg})))
	use_quatern=1;
      end
      if ~isempty(findstr('qform44',varargin{currarg}))
	use_qform=1;
      end
      if strcmp(varargin{currarg},'datatype')
	hdr.matlab_datatype=cbiNiftiDatatype2Matlab(hdr.datatype);
	hdr.bitpix=8*cbiSizeofNifti(hdr.matlab_datatype);
      end
      if strcmp(varargin{currarg},'matlab_datatype')
	hdr.datatype=cbiMatlabDatatype2Nifti(hdr.matlab_datatype);
	hdr.bitpix=8*cbiSizeofNifti(hdr.matlab_datatype);
      end
    end
  end
  if (use_qform)
    use_quatern=0;
  end
  if (use_sform)
    use_srow=0;
  end    

  if (~use_qform & ~use_quatern)
    use_qform=1;
  end
  if (~use_sform & ~use_srow)
    use_sform=1;
  end 
  
  % Set data size, if any
  if (input_data)
    s=size(varargin{input_data});    
    l=length(s);
    hdr.dim(:)=0;
    hdr.dim(1)=l;
    hdr.dim(2:l+1)=s;
    if (l>2)
      hdr.slice_end=s(3)-1; % dim3 (Z) is slice direction
    else
      hdr.slice_end=0;
    end
    % Set pixdim fields, if not set
    for n=2:l+1
      if (hdr.pixdim(n)==0)
	hdr.pixdim(n)=1;
      end
    end
  end
    
  if (use_quatern) 
    % Make sure quaternions are valid
    quatern=[hdr.quatern_b hdr.quatern_c hdr.quatern_d];
    pixdim=hdr.pixdims(2:4);
    qfac=hdr.pixdims(1);
    qoffs=[hdr.qoffset_x;hdr.qoffset_y;hdr.qoffset_z];
    if (length(quatern)==3 & length(qoffs)==3 & prod(pixdim)~=0)
      % Calculate qform44
      hdr.qform44=cbiQuaternionToHomogeneous( quatern, pixdim, qfac, qoffset );
    else
      % Else set all to null values
      hdr.quatern_b=0;
      hdr.quatern_c=0;
      hdr.quatern_d=0;
      hdr.qoffset_x=0;
      hdr.qoffset_y=0;
      hdr.qoffset_z=0;
      hdr.qform44=[];
    end
  end
  if (use_qform)
    % Make sure qform44 is valid
    if (~isempty(hdr.qform44) & ~((hdr.qform_code==0) & isequal(hdr.qform44,eye(4))) & ~isequal(hdr.qform44,zeros(4,4)) & size(hdr.qform44)==[4,4])
      hdr=cbiSetNiftiQform( hdr, hdr.qform44 );
    else
      % Else set all to null values
      hdr.quatern_b=0;
      hdr.quatern_c=0;
      hdr.quatern_d=0;
      hdr.qoffset_x=0;
      hdr.qoffset_y=0;
      hdr.qoffset_z=0;
      hdr.qform44=[];
    end
  end
  if (use_srow)
    % Make sure srows are valid
    if (length(hdr.srow_x)==4 & length(hdr.srow_y)==4 & length(hdr.srow_z)==4 )      
      % Calculate sform44
      hdr.sform44=eye(4);
      hdr.sform44(1,:)=hdr.srow_x;
      hdr.sform44(2,:)=hdr.srow_y;
      hdr.sform44(3,:)=hdr.srow_z;      
    else
      % Else set all to null values
      hdr.srow_x=zeros(1,4);
      hdr.srow_y=zeros(1,4);
      hdr.srow_z=zeros(1,4);
      hdr.sform44=[];
    end
  end
  if (use_sform)
    % Make sure sform44 are valid
    if (~isempty(hdr.sform44) & ~((hdr.sform_code==0) & isequal(hdr.sform44,eye(4))) & ~isequal(hdr.sform44,zeros(4)) & size(hdr.sform44)==[4,4])
      % Calculate srows
      hdr.srow_x=hdr.sform44(1,:);
      hdr.srow_y=hdr.sform44(2,:);
      hdr.srow_z=hdr.sform44(3,:);
    else
      % Else set all to null values
      hdr.srow_x=zeros(1,4);
      hdr.srow_y=zeros(1,4);
      hdr.srow_z=zeros(1,4);
      hdr.sform44=[];
    end
  end
  
  
  % Set qform_code, sform_code according to header
  % qform44: 
  if (~isempty(hdr.qform44))
    hdr.qform_code=1;
  else
    hdr.qform_code=0;
  end
  % sform44:
  if (~isempty(hdr.sform44))
    % only change sform_code if it hasn't been set
    % otherwise leave it unchanged (since it could be =3,not =1)
    if isequal(hdr.sform_code,0) || isempty(hdr.sform_code)
      hdr.qform_code=1;
    end
  else
    hdr.sform_code=0;
  end

  % file name: set magic number
  if (isfield(hdr,'hdr_name'))
    [pathstr,bname,ext]=fileparts(hdr.hdr_name);
    switch (ext)
     case '.nii'  
       hdr.single_file=1;
       hdr.magic='n+1';
     case {'.hdr','.img'}
      hdr.single_file=0;
      hdr.magic='ni1';      
     otherwise
      disp('Warning: unrecognized file name field.');
    end
  end
  
  % Truncate char arrays
  if (length(hdr.aux_file)>24)
    hdr.aux_file=hdr.aux_file(1:24);
  else
    hdr.aux_file=[hdr.aux_file char(zeros(1,24-length(hdr.aux_file)))];
  end
  if (length(hdr.descrip)>80)
    hdr.descrip=hdr.descrip(1:80);
  else
    hdr.descrip=[hdr.descrip char(zeros(1,80-length(hdr.descrip)))];
  end
  if (length(hdr.intent_name)>16)
    hdr.intent_name=hdr.intent_name(1:16);
  else
    hdr.intent_name=[hdr.intent_name char(zeros(1,16-length(hdr.intent_name)))];
  end
  
  return
  
  
function h=generateEmptyNiftiHeader(h)
  if (~exist('h'))
    h=[];
  end
  % For every header field, check if it exists, otherwise add it
  h.sizeof_hdr=348;
  if (~isfield(h,'dim_info'))
    h.dim_info=' ';
  end
  if (~isfield(h,'dim'))
    h.dim=zeros(8,1);
  end
  if (~isfield(h,'intent_ps'))
    h.intent_ps=zeros(3,1);
  end    
  if (~isfield(h,'intent_code'))
    h.intent_code=0;
  end
  if (~isfield(h,'endian'))
    h.endian='native';
  end
  if (~isfield(h,'datatype'))
    h.datatype=16;
  end
  h.matlab_datatype=cbiNiftiDatatype2Matlab(h.datatype);
  h.bitpix=8*cbiSizeofNifti(h.matlab_datatype);
  if (~isfield(h,'slice_start'))
    h.slice_start=0;
  end
  if (~isfield(h,'slice_end'))
    h.slice_end=0;
  end
  if (~isfield(h,'pixdim'))
    h.pixdim=zeros(8,1);
  end
  
  if (~isfield(h,'vox_offset'))
    h.vox_offset=0;
  end
  if (~isfield(h,'scl_slope'))
    h.scl_slope=1;
  end
  if (~isfield(h,'scl_inter'))
    h.scl_inter=0;
  end
  
  if (~isfield(h,'slice_code'))
    h.slice_code=0;
  end
  
  if (~isfield(h,'xyzt_units'))
    % Note, this used to be wrong before 6/17/2015 and
    % was set to 18 which is not secs but ms. Fixed jg
    h.xyzt_units=10; % NIFTI_UNITS_MM | NIFTI_UNITS_SEC
  end
  if (~isfield(h,'cal_max'))
    h.cal_max=0;
  end
  if (~isfield(h,'cal_min'))
    h.cal_min=0;
  end
  if (~isfield(h,'slice_duration'))
    h.slice_duration=0;
  end
  if (~isfield(h,'toffset'))
    h.toffset=0;
  end
  if (~isfield(h,'descrip'))
    h.descrip=char(zeros(1,80));
  end
  if (~isfield(h,'aux_file'))
    h.aux_file=char(zeros(1,24));
  end
  if (~isfield(h,'qform_code'))
    h.qform_code=0;
  end
  if (~isfield(h,'sform_code'))
    h.sform_code=0;
  end
  if (~isfield(h,'quatern_b'))
    h.quatern_b=0;
  end
  if (~isfield(h,'quatern_c'))
    h.quatern_c=0;
  end
  if (~isfield(h,'quatern_d'))
    h.quatern_d=0;
  end
  if (~isfield(h,'qoffset_x'))
    h.qoffset_x=0;
  end
  if (~isfield(h,'qoffset_y'))
    h.qoffset_y=0;
  end
  if (~isfield(h,'qoffset_z'))
    h.qoffset_z=0;
  end
  if (~isfield(h,'srow_x'))
    h.srow_x=[0 0 0 0];
  end
  if (~isfield(h,'srow_y'))
    h.srow_y=[0 0 0 0];
  end
  if (~isfield(h,'srow_z'))
    h.srow_z=[0 0 0 0];
  end
  if (~isfield(h,'intent_name'))
    h.intent_name=char(zeros(1,16));
  end
  if (~isfield(h,'magic'))
    h.magic='';
  end
  if (~isfield(h,'qform44'))
    h.qform44=[];
  end
  if (~isfield(h,'sform44'))
    h.sform44=[];
  end    
  if (~isfield(h,'intent_ps'))
    h.intent_ps=zeros(3,1);
  end
  return
  
      
  
  
  
