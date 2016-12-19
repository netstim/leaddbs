function [res, errStr, res2]= maskstruct_query(varargin)
%function [res, errStr, res2]= maskstruct_query(maskStruct, commandStr[, op1[, op2[... [, opN]]]])
%
%   General method query a infomations and data from a maskStruct
%
%   maskStruct: should contain a valid maskStruct Ver2
%   command:    {}
%   op[1..end] the parameter specified by the command string
%
% 'sizeAy'
% 'maskNo'
% 'maskNames'
% 'getMaskId' nameStr (string)
% 'getMaskIdx' nameStr (string|int), [extendFlag]
% 'getMask' nameStr (string|int), [extendFlag]
% 'getMRMask'   nameStr (string|int), [extendFlag]
% 'testMRDat'   data (mrStruct|dtdStruct)
%           res: -1 - not combinable (different size and no edges)
%                 0 - different size but coregisteration is possible
%                 1 - all fine; same size and same orientation
%           res2: true  - edges information was available in maskstruct
%                 false - no valid edges information is available in maskstruct
% 'getMassPoint nameStr (string|int), [nearstFlag] ('yes')
% 
% 
% return values:
%   res: The result of the method or if an error occured res is empty
%   errStr: a message string, specifying the error which occured
%   res2: additional output argument, depending on the command
% 
% Description of the different commands: 
% Information about the maskStruct: 
%  'fiberNames'returns the names (identifier) of all fiber subsets (fiberStruct) in th ftrStruct
%              res: cell string containing the fiber subset names
% 
%
%
% Bjoern W. Kreher
% 12/07
%
% UNIX


res= [];    res2= [];    errStr= '';
%% check input param
if length(varargin) < 2
    errStr= sprintf('%s(error): There have to be at least two parameters', mfilename);
    return;
end

if (~isempty(varargin)) && maskstruct_istype(varargin{1})
    maskStruct= varargin{1};
else
    errStr= sprintf('%s(error): First param have to be maskStruct Ver2', mfilename);
    return;
end

if ~ischar(varargin{2})
    errStr= sprintf('%s(error): Command have to be a string', mfilename);
    return;
else
    commandStr= varargin{2};
end

%% move input param
paramMax= 10;
param= cell(paramMax, 1);
for i= 1:paramMax
    if length(varargin) < (i + 2)
        param{i}= [];
    else
        param{i}= varargin{i + 2};
    end
end


%% begin of command switching part
if strcmp(commandStr, 'sizeAy')
    [res, errStr]= local_getVolSize(maskStruct);
elseif strcmp(commandStr, 'maskNo')
    [res, errStr]= local_maskNo(maskStruct);
elseif strcmp(commandStr, 'maskNames')
    [res, errStr]= local_getMaskNames(maskStruct);
elseif strcmp(commandStr, 'isMaskName')
    if (ischar(param{1}) || isscalar(param{1})) 
        nameStr= param{1};  
    else
        errStr= sprintf('%s(isMaskName): invalid data', mfilename);
        return
    end    
    [res, errStr]= local_isMaskName(maskStruct, nameStr);
elseif strcmp(commandStr, 'getMaskIdx')
    if (ischar(param{1}) || isscalar(param{1})) && (isempty(param{2}) || ischar(param{2}))
        nameStr= param{1};  extendFlag= param{2}; %'noExtend'
    else
        errStr= sprintf('%s(getMaskIdx): invalid data', mfilename);
        return
    end    
    [res, errStr]= local_getMaskIdx(maskStruct, nameStr, extendFlag);
elseif strcmp(commandStr, 'getMaskVc')
    if (ischar(param{1}) || isscalar(param{1})) 
        nameStr= param{1};  
    else
        errStr= sprintf('%s(getMaskVc): invalid data', mfilename);
        return
    end    
    [res, errStr]= local_getMaskVc(maskStruct, nameStr);
elseif strcmp(commandStr, 'getMaskName')
    if isscalar(param{1}) && (param{1} == round(param{1}))
        id= param{1};
    else
        errStr= sprintf('%s(getMaskName): invalid data', mfilename);
        return
    end    
    [res, errStr]= local_getMaskName(maskStruct, id);
elseif strcmp(commandStr, 'getMaskId')
    if ischar(param{1}) 
        nameStr= param{1};
    else
        errStr= sprintf('%s(getMaskId): invalid data', mfilename);
        return
    end    
    [res, errStr]= local_getMaskId(maskStruct, nameStr);
elseif strcmp(commandStr, 'getMask')
    if (ischar(param{1}) || isscalar(param{1})) && (isempty(param{2}) || ischar(param{2}))
        nameStr= param{1};  extendFlag= param{2}; %'noExtend'
    else
        errStr= sprintf('%s(getMask): invalid data', mfilename);
        return
    end    
    [res, errStr]= local_getMask(maskStruct, nameStr, extendFlag);
elseif strcmp(commandStr, 'getMRMask')
    if (ischar(param{1}) || isscalar(param{1})) && (isempty(param{2}) || ischar(param{2}))
        nameStr= param{1};  extendFlag= param{2}; %'noExtend'
    else
        errStr= sprintf('%s(getMRMask): invalid data', mfilename);
        return
    end    
    [res, errStr]= local_getMRMask(maskStruct, nameStr, extendFlag);
elseif strcmp(commandStr, 'testMRDat')
    if dtdstruct_istype(param{1}) 
        dtdStruct= param{1};
    else
        errStr= sprintf('%s(testMRDat): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_testMRDat(maskStruct, dtdStruct);
elseif strcmp(commandStr, 'getMassPoint')
    if (ischar(param{1}) || isscalar(param{1})) && (isempty(param{2}) || ischar(param{2}))
        nameStr= param{1};  extendFlag= param{2}; %'yes'
    else
        errStr= sprintf('%s(getMassPoint): invalid data', mfilename);
        return
    end    
    [res, errStr]= local_getMassPoint(maskStruct, nameStr, extendFlag);
else
    errStr= sprintf('%s(varagrin): command ''%s'' is not implemented', mfilename, commandStr);
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, orInfo]= local_testMRDat(maskStruct, dtdStruct)
errStr= '';
res= -1;
orInfo= false;


orInfo= isequal(size(mrstruct_query(maskStruct.mrsProp, 'edges')), [4 4]);
[res, errStr]= dtdstruct_query(dtdStruct, 'cmpDimensions', ...
    mrstruct_init('volume', false(maskStruct.sizeAy), maskStruct.mrsProp));



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_isMaskName(maskStruct, nameStr)
errStr= '';
res= any(strcmp(maskStruct.maskNamesCell, nameStr));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getMaskNames(maskStruct)
errStr= '';
res= maskStruct.maskNamesCell;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getVolSize(maskStruct)
errStr= '';
res= maskStruct.sizeAy;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_maskNo(mStruct)
errStr= '';
res= length(mStruct.maskCell);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getMaskName(mStruct, id)
res= []; errStr= '';

if (id <= 0) || (id > length(mStruct.maskNamesCell))
    errStr= 'maskstruct_query(local_getMaskName): invalid mask id';
else
    res= mStruct.maskNamesCell{id};
end
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getMaskId(mStruct, nameStr)
res= []; errStr= '';

[res, errStr]= private_name2id(mStruct, nameStr);
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getMaskIdx(mStruct, nameStr, extFlag)
res= []; errStr= '';

[id, errStr]= private_name2id(mStruct, nameStr);
if isempty(id)
    return;
end

if ~isempty(mStruct.maskCell{id})
    res= find(mStruct.maskCell{id});
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getMassPoint(mStruct, nameStr, extFlag)
res= []; errStr= '';

[id, errStr]= private_name2id(mStruct, nameStr);
if isempty(id)
    return;
end

idx= find(mStruct.maskCell{id});
if isempty(idx)
    errStr= 'maskstruct_query::local_getMassPoint(error): Mask is empty';
    return;
end

dVc= reshape_index(idx, mStruct.sizeAy);
res= mean(dVc, 1);

if strcmp(extFlag, 'yes')
    [val, idx]= min((dVc(:, 1) - res(1)).^2 + (dVc(:, 2) - res(2)).^2 + (dVc(:, 3) - res(3)).^2);
    res= dVc(idx, :);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getMaskVc(mStruct, nameStr)
res= []; errStr= '';

[id, errStr]= private_name2id(mStruct, nameStr);
if isempty(id)
    return;
end

idx= find(mStruct.maskCell{id});
if isempty(idx)
    errStr= 'maskstruct_query::local_getMaskVc(error): Mask is empty';
    return;
end

res= reshape_index(idx, mStruct.sizeAy);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getMask(mStruct, nameStr, extFlag)
res= []; errStr= '';

[id, errStr]= private_name2id(mStruct, nameStr);
if isempty(id)
    return;
end

res= mStruct.maskCell{id};

if isempty(res) && ~strcmp(extFlag, 'noExtend')
    if ~isempty(mStruct.sizeAy)
        res= false(mStruct.sizeAy);
    else
        errStr= 'maskstruct_query(local_getMask): can not extend data. No volume information';
        return;
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getMRMask(mStruct, nameStr, extFlag)
res= []; errStr= '';

[dat, errStr]= local_getMask(mStruct, nameStr, extFlag);
if isempty(errStr)
   res= mrstruct_init('volume', double(dat), mStruct.mrsProp);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  START:
%       [id, errStr]= private_name2id(mStruct, nameStr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [id, errStr]= private_name2id(mStruct, nameStr)
id = []; 
errStr= '';

if isnumeric(nameStr)
    id= nameStr;
    if (id < 1) || (id > length(mStruct.maskNamesCell)) && (round(id) ~= id)
        id= [];
        errStr= 'maskstruct_modify(private_name2id): invalid mask id';
    end
    return
end

if ~ischar(nameStr)
    id= [];
    errStr= 'maskstruct_modify(private_name2id): invalid mask id';
    return;
end

id= find(strcmp(mStruct.maskNamesCell, nameStr));

if isempty(id)
    errStr= 'maskstruct_modify(private_name2id): invalid mask name';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  END: private_name2id
