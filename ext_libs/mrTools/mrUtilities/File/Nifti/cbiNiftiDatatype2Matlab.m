function matlab_datatype=cbiNiftiDatatype2Matlab(nifti_datatype);
% function matlab_datatype=cbiNiftiDatatype2Matlab(nifti_datatype);
%
% Translates the numeric (integer) datatype code of a Nifti header to the corresponding Matlab datatype
  
switch(nifti_datatype)
 % ANALYZE & NIFTI 
 case 1
  matlab_datatype='binary';
 case 2
  matlab_datatype='uint8';
 case 4 
  matlab_datatype='int16';
 case 8
  matlab_datatype='int32';  
 case 16
  matlab_datatype='float32';
 case 32
  matlab_datatype='complex64';
 case 64
  matlab_datatype='double';
 case 128
  matlab_datatype='RGB';
  % NIFTI specific (long double and long double pair only properly supported on 128-bit systems)
 case 256
  matlab_datatype='int8';
 case 512
  matlab_datatype='uint16';
 case 768
  matlab_datatype='uint32';
 case 1024
  matlab_datatype='int64';
 case 1280
  matlab_datatype='uint64';
 case 1536
  matlab_datatype='float128';  
 case 1792
  matlab_datatype='complex128';  
 case 2048
  matlab_datatype='complex256';  
end

