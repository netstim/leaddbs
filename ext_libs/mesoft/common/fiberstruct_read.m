function [res, errStr, fName]= fiberstruct_read(fName)
%
%function [res, errStr, fName]= fiberstruct_read(fName)
%
%
%  Bjoern W. Kreher
%  08/02
%
%  UNIX
%

errStr= '';

if nargin == 0 
    [fileStr,dirStr]=uigetfile('*.mat','Load a fiberStruct');
    fName= strcat(dirStr, '/', fileStr);
end
res= open(fName);
