function [res, errStr]= ftrstruct_init(typeStr)
%   function [res, errStr]= ftrstruct_init(typeStr)
%
%  Creates a valid and empty ftrStruct or fiberStruct
%
%   typeStr: {'ftrStruct' | 'fiberStruct'}
%           default is 'ftrStruct'
%
%   return values
%     res:  the empty ftrStruct or fiberStruct
%     errStr: If an error occured, errStr is identify the error
%
%
% Bjoern W. Kreher
% 11/02
%
% UNIX

%       typeStr= {'DTD' | 'multiDTD' | 'fiber'};


errStr= '';
res= [];
if exist('typeStr')
    res= local_getBaseStruct;
    if strcmp(typeStr, 'DTD')
        res.dtdType= typeStr;
    elseif strcmp(typeStr, 'multiDTD')
        res.dtdType= typeStr;
    elseif (strcmp(typeStr, 'fiber') || strcmp(typeStr, 'fiberStruct'))
        res= local_getFiberStruct;
    elseif strcmp(typeStr, 'ftrStruct')
    else
        errStr= 'ftrstruict_init: no valid typeStr';
        res= [];
    end
else
    res= local_getBaseStruct;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%
%  START:
%       [res, errStr]= local_getBaseStruct
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getBaseStruct

res.connectCell=    {};
res.curveSegCell=   {};
res.posSegCell=     {};

res.fiber=          {};

res.dtdType=        '';
res.vox=            [];
res.patient=        '';
res.algoName=       '';
res.trackDate=      '';
res.trackParam=     [];
res.logData=        [];
res.hMatrix=        [];
res.user=           [];
res.version=        'V1.1';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getBaseStruct
%

%
%
%  START:
%       [res, errStr]= local_getFiberStruct
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getFiberStruct

res.name=           '';
res.curveID=        [];
res.roiName=        {};
res.user=           [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getFiberStruct
%


