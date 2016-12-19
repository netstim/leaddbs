function [res, errStr, verStr]= probstruct_istype(dataIn)
%Tests if the dataIn is a valid probStruct and extract the VersionTag of the probStruct
%   function [res, errStr, verStr]= probstruct_istype(dataIn)
%
%   test function returns 1 if dataIn is a valid probStruct 
%
%   dataIn:    arbitary data .. 
%   fiberType: {'Ver1' | 'Ver2'} the default value is 'Ver2'
%
%   return values
%   res:    true if dataIn is an probStruct of Version 'V2' else it wil be false
%   errStr: If an error occured, errStr is identify the error
%   verStr: Returns the versionTag 'V2'; 'V1' or '' 
%
% Bjoern W. Kreher
% 12/07
%
% UNIX

res= false; errStr= ''; verStr= ''; 

%% Test input parameter
if nargin < 1
    errStr= sprintf('%s(error): function needs at list one parameter', mfilename);
    return;
end

if ~mrstruct_istype(dataIn)
    errStr= sprintf('%s(error): data is not mrStruct', mfilename);
    return;
end



if local_test_probStruct_V1(dataIn.user)
    verStr= 'V1';
else
    [res, errStr]= local_test_probStruct_V2(dataIn.user);
    if res
        verStr= 'V2';
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_test_probStruct_V1(user)
res= false; errStr= '';

if ~all(isfield(user, {'mapName', 'algName', 'genDate', 'versionLog', 'operationLog'}))
    errStr= 'missing field in user';
    return
end

if ~(ischar(user.mapName) && ischar(user.algName) && ischar(user.versionLog) && ischar(user.genDate) && isstruct(user.operationLog))
    errStr= 'wrong type in user';
    return
end

if ~isfield(user.operationLog, 'operationLog') || ~ischar(user.operationLog.operationLog)
    errStr= 'operationLog in operation Log';
    return
end

res= true;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_test_probStruct_V2(user)
res= false; errStr= '';

if ~isfield(user, 'probInfo')
    errStr= 'no probInfo field in user';
    return
end

if ~all(isfield(user.probInfo, {'mapName', 'algType', 'dataType', 'operationLog', 'histStruc', 'logHistStruc', 'MinMaxStruc'}))
    errStr= 'missing field in probInfo';
    return
end

if ~all(isfield(user.probInfo.operationLog, {'type', 'date', 'params'}))
    errStr= 'missing field in operationLog';
    return
end

res= true;



