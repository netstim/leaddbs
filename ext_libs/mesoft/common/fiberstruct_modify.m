function [res, errStr]= fiberstruct_modify(varargin)
%
%   function [res, errStr]= fiberstruct_query(fiberStruct, command, parm1[, parm2[, ...[, parmN]]])
%
%   command:    ['trans' | ]
%           setOM: parm1= hmMatrix
%
%  Bjoern W. Kreher
%  08/02
%
%  UNIX
%

res= [];    errStr= '';

if ischar(varargin{1}) && exist(varargin{1}, 2)
    fiberStruct= fiberstruct_read(varargin{1});
elseif fiberstruct_istype(varargin{1})
    fiberStruct= varargin{1};
else
    errStr= 'Error in function [res, errStr]= fiberstruct_modify(varagrin): first param have to be fiberStruct';
    return;
end

if ~ischar(varargin{2})
    errStr= 'Error in function [res, errStr]= fiberstruct_modify(varagrin): command have to be a string';
    return;
end
commandStr= varargin{2};

if strcmp(commandStr, 'trans')
    orM= varargin{3};
    [res, errStr]= local_transCurve(fiberStruct, orM);
else
    errStr= sprintf('Error in function [res, errStr]= fiberstruct_modify(varagrin): command ''%s'' is not implemented', commandStr);
end






function [res, errStr]= local_transCurve(fiberStruct, orM)

res.curvesData= {}; errStr= '';

for i= 1:length(fiberStruct.curvesData)
    res.curvesData{i}.roiVertex= fiberStruct.curvesData{i}.roiVertex;
    res.curvesData{i}.curves= {};
    res.curvesData{i}.pos= {};
    res.curvesData{i}.error= {};
    for k= 1:length(fiberStruct.curvesData{i}.curves)
        res.curvesData{i}.curves{k}= hm_mult(orM, fiberStruct.curvesData{i}.curves{k}-.5)+.5;
        res.curvesData{i}.pos{k}= hm_mult(orM, fiberStruct.curvesData{i}.pos{k}-.5)+.5;
        res.curvesData{i}.error{k}= hm_mult(orM, fiberStruct.curvesData{i}.error{k});
    end
end


