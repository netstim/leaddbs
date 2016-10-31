function nifti_datatype=cbiMatlabDatatype2Nifti(matlab_datatype);
% function nifti_datatype=cbiMatlabDatatype2Nifti(matlab_datatype);
%
% Translates the Matlab datatype into a numeric (integer) Nifti datatype code
  
switch(matlab_datatype)
 % ANALYZE & NIFTI 
 case 'binary'
  nifti_datatype=1;
 case 'uint8'
  nifti_datatype=2;
 case {'int16','short'}
  nifti_datatype=4;
 case {'int32','int'}
  nifti_datatype=8;
 case 'float32'
  nifti_datatype=16;
 case 'complex64'
  nifti_datatype=32;
 case 'double';
  nifti_datatype=64;
 case 'RGB'
  nifti_datatype=128;
  % NIFTI specific (long double and long double pair only properly supported on 128-bit systems)
 case 'int8'
  nifti_datatype=256;
 case {'uint16','ushort'}
  nifti_datatype=512;
 case {'uint32','uint'}
  nifti_datatype=768;
 case 'int64'
  nifti_datatype=1024;
 case 'uint64'
  nifti_datatype=1280;
 case 'float128'  
  nifti_datatype=1536;
 case 'complex128'  
  nifti_datatype=1792;
 case 'complex256'
  nifti_datatype=2048;
end

