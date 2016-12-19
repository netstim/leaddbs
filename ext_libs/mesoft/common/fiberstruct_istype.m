function [res, errStr]= fiberstruct_istype(fiberStruct)
%
%  function [res, errStr]= fiberstruct_istype(fiberStruct)
%
%  Bjoern W. Kreher
%  08/02
%
%  UNIX
%

res= 0;
if isstruct(fiberStruct)
    str= fieldnames(fiberStruct);
    if length(str) == 1
        res= strcmp(str,'curvesData');
    end
end

