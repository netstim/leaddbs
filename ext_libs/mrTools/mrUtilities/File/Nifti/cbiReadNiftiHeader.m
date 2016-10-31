function hdr = cbiReadNiftiHeader(fname)
% hdr = cbiReadNiftiHeader(fname)
% 
% Loads the header from a NIFTI-1 file.  
% The header struct contains all stored fields plus 
% some additional extracted fields. These are:
% .hdr_name    - name of header file
% .img_name    - name of image file
% .is_analyze  - 1 if valid Analyze, but not a valid NIFTI file, 0 otherwise
% .single_file - flag for single file (1 for .nii, 0 for .hdr/.img)
% .endian      - endianness of file
% .matlab_datatype - data type in Matlab string format (see fread)
% .qform44    - homogeneous transform matrix extracted from qform data (rigid body - native coordinate system) ([] for Analyze) 
% .sform44    - homogeneous transform matrix extracted from sform data (affine - optional alignment transform) ([] for Analyze) 
%
% fname can have any legal extension (.hdr, .img, .nii). No support for compressed data - will bail out.

% Check name and file type.
[pathstr,bname,ext]=fileparts(fname);

switch (ext)
 case '.nii'  
  hdr.single_file=1;
  hdr.hdr_name=fname;
  hdr.img_name=fname;
 case {'.hdr','.img'}
  hdr.single_file=0;
  hdr.hdr_name=fullfile(pathstr,[bname '.hdr']);
  hdr.img_name=fullfile(pathstr,[bname '.img']);
 case '.gz','.Z' % zipped
  disp(sprintf('(cbiReadNifitHeader) No support for zipped NIFTI-1 format under Matlab.'));
  hdr = [];
  return
 otherwise
  disp(sprintf('(cbiReadNiftiHeader) %s is not a valid NIFTI-1 file name extension. Legal values are .nii, .hdr, .img',fname));
  hdr = [];
  return;
end

% Find the endian-ness of the file
hdr.endian='b'; % open file in big-endian
fid=fopen(hdr.hdr_name,'r',hdr.endian);
if fid==-1,hdr=[];return,end
% check if this gives the correct header size - if not use little-endian
testval = fread(fid,1,'int32');
if ~isequal(testval,348)
  fclose(fid);
  hdr.endian='l';
  fid=fopen(hdr.hdr_name,'r',hdr.endian);
  testval = fread(fid,1,'int32');
  if ~isequal(testval,348)
    disp(sprintf('(cbiReadNiftiHeader) Incorrect header size (should be 348 bytes) for file %s',fname));
    hdr = [];
    return
  end
end

%% --- was header_key substruct ---
hdr.sizeof_hdr = testval;
dummy = fread(fid,35,'char');	% ditch the remaining initial header stuff
hdr.dim_info = fread(fid,1,'char');
%% --- was image_dimension substruct ---
hdr.dim = fread(fid,8,'int16');
hdr.intent_ps = fread(fid,3,'float');
hdr.intent_code = fread(fid,1,'int16');
hdr.datatype = fread(fid,1,'int16');
hdr.bitpix = fread(fid,1,'int16');
hdr.slice_start = fread(fid,1,'int16');
hdr.pixdim = fread(fid,8,'float');
hdr.vox_offset = fread(fid,1,'float');
hdr.scl_slope = fread(fid,1,'float');
hdr.scl_inter = fread(fid,1,'float');
hdr.slice_end = fread(fid,1,'int16');
hdr.slice_code = fread(fid,1,'char');
hdr.xyzt_units = fread(fid,1,'char');
hdr.cal_max = fread(fid,1,'float');
hdr.cal_min = fread(fid,1,'float');
hdr.slice_duration = fread(fid,1,'float');
hdr.toffset = fread(fid,1,'float');
dummy = fread(fid,2,'int32');
%% --- was data_history substruct ---
hdr.descrip = char(transpose(fread(fid,80,'char')));
hdr.aux_file = char(transpose(fread(fid,24,'char')));

hdr.qform_code = fread(fid,1,'int16');
hdr.sform_code = fread(fid,1,'int16');
hdr.quatern_b = fread(fid,1,'float');
hdr.quatern_c = fread(fid,1,'float');
hdr.quatern_d = fread(fid,1,'float');
hdr.qoffset_x = fread(fid,1,'float');
hdr.qoffset_y = fread(fid,1,'float');
hdr.qoffset_z = fread(fid,1,'float');
hdr.srow_x    = fread(fid,4,'float')';
hdr.srow_y    = fread(fid,4,'float')';
hdr.srow_z    = fread(fid,4,'float')';
hdr.intent_name  = char(transpose(fread(fid,16,'char')));
hdr.magic        = char(transpose(fread(fid,4,'char')));

fclose(fid);

% Set data type
hdr.matlab_datatype=cbiNiftiDatatype2Matlab(hdr.datatype);

hdr.is_analyze=(~hdr.single_file & ~strfind(hdr.magic,'ni1'));

% Extract xform and sform
if (hdr.is_analyze)
  hdr.qform44=[];
  hdr.sform44=[];
else
  hdr.qform44=cbiQuaternionToHomogeneous(hdr);
  hdr.sform44=eye(4);
  if (hdr.sform_code>0)
    hdr.sform44(1,:)=hdr.srow_x;
    hdr.sform44(2,:)=hdr.srow_y;
    hdr.sform44(3,:)=hdr.srow_z;
  end
end
