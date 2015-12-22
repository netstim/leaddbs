function [res, errStr]= maskstruct_istype(dataIn, maskType)
%function [res, errStr]= maskstruct_istype(dataIn, maskType)
%
%   test function returns 1 if dataIn is a valid maskStruct of 
%   the version specified in maskType. 
%
%   dataIn:    arbitary data .. may be a ftrStruct or a fiberStruct
%   fiberType: {'Ver1' | 'Ver2'} the default value is 'Ver2'
%
%   return values
%   res:    1 if dataIn is an ftrStruct 0 else
%   errStr: If an error occured, errStr is identify the error
%
% Bjoern W. Kreher
% 12/07
%
% UNIX

res= false; errStr= ''; 

%% Test input parameter
if ~exist('dataIn')
    errStr= sprintf('%s(error): function needs at list one parameter', mfilename);
    return;
end
if ~exist('maskType') || isempty(maskType)
    maskType= 'Ver2';
end
if ~ischar(maskType)
    errStr= sprintf('%s(error): parameter maskType must be of type string', mfilename);
    return;
end

%% split between different versions
if strcmp(maskType, 'Ver1')
    [res, errStr]= local_istypeV1(dataIn);
elseif strcmp(maskType, 'Ver2')
    [res, errStr]= local_istypeV2(dataIn);
else
    errStr= sprintf('%s(error): The Version ''%s'' is not supported yet', mfilename, maskType);
    return;    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_istypeV1(dataIn)
res= false; errStr= '';

if ~isstruct(dataIn)
    return;
end
dataCell= struct2cell(dataIn);
sizeOrgAy= [];
for i= 1:length(dataCell)
    if ~isempty(dataCell{i})
       if mrstruct_istype(dataCell{i})
           if isempty(sizeOrgAy)
               sizeOrgAy= mrstruct_query(dataCell{i}, 'sizeAy'); %% Setze Gr�sse
           elseif ~isequal(sizeOrgAy, mrstruct_query(dataCell{i}, 'sizeAy'))
               return %% falsche Gr�sse
           end
       else 
           return %% kein mrStruct
       end
    end
end
res= true;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_istypeV2(dataIn)
res= false; errStr= '';

if ~isstruct(dataIn)
    return
end

if ~all(isfield(dataIn, {'maskCell'; 'maskNamesCell'; 'sizeAy'; 'mrsProp'; 'user'; 'version'}))
    return;
end
if length(dataIn.maskCell) ~= length(dataIn.maskNamesCell)
    return;
end

res= true;



