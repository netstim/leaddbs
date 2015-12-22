function [res, errStr]= maskstruct_init(mrStruct)
%   function [res, errStr]= maskstruct_init(mrStruct)
%
%  Creates a valid and empty maskStruct
%
%   return values
%     res:  the empty maskStruct
%     errStr: If an error occured, errStr is identify the error
%
%
% Bjoern W. Kreher
% 12/07
%
% UNIX


res= [];   errStr= '';

if exist('mrStruct','var') && ~isempty(mrStruct)
    if dtdstruct_istype(mrStruct)
        res= maskstruct_modify(local_getBaseStruct, 'mrStructProb', dtdstruct_query(mrStruct, 'mrStructProb'));
        res= maskstruct_modify(res, 'setVolumeSize', dtdstruct_query(mrStruct, 'sizeAy'));        
    else
        errStr= 'maskstruct_init(error): Arg. is not of type mr- or dtdStruct';
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
errStr= '';
res.maskCell=       {};
res.maskNamesCell=  {};
res.sizeAy=         [];

res.mrsProp=        [];
res.user=           [];
res.version=        'V2.0';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getBaseStruct
%


