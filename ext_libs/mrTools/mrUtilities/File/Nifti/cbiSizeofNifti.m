function varargout=cbiSizeofNifti(datatype);
% [nbytes,possnbytes]=cbiSizeofNifti(datatype);
% 
% Returns the number of bytes for a given Nifti data type
% in string format (see fread for details).
% Returns [] for unrecognized data types.
%
% WARNING: some data types are not unambiguously defined.
% For these, the alternate number of bytes is output in 
% possnbytes, if 2 output arguments are given.
% These data types are:
% Type     nbytes possnbytes
% long     4      8
% ulong    4      8

switch (datatype)
 case {'uchar','schar','char','uint8','int8'}
  nbytes=1;
  possnbytes=nbytes;
 case {'int16','uint16','short','ushort'}
  nbytes=2;
  possnbytes=nbytes;
 case {'int32','uint32','single','float32','int','uint'}
  nbytes=4;
  possnbytes=nbytes;
 case {'int64','uint64','float64','double','complex64'}
  nbytes=8;
  possnbytes=nbytes;
 case {'long','ulong'}
  nbytes=4;
  possnbytes=8;
 case {'float128','complex128'}
  nbytes=16;
  possnbytes=nbytes;
 case 'complex256'
  nbytes=32;
  possnbytes=nbytes;
 otherwise
  nbytes=[];
  possnbytes=[];
end
varargout{1}=nbytes;
if (nargout==2)
  varargout{2}=possnbytes;
end

