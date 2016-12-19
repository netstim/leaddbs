function [res, errStr]= probstruct_init(mrStruct, mapName, dataType, algType, opLog)
%creates a new probStruct object
%  [res, errStr]= probstruct_init(mrStruct, mapName, dataType, algType, opLog)
%
%  Creates a valid probStruct V2 from a mrStruct 
%       
%       mrStruct - (mrStruct) a containing also the map data 
%       mapName  - (String) A name specifying the map
%       dataType - (String) <'2ModMult' | 'Mult' | '2ModPMap' | 'PMap'>
%                   specifying the type of map
%                   'PMap' containing a measure of connectivity to a seed
%                       point or similar
%                   'PMult' result of a multiplication between two ('PMap'*'PMap' -> 'PMult')
%                   '2ModPMap' containing a two modal measure of connectivity to a seed
%                       point or similar (eg. comming form probRandExt)
%                   '2ModPMult' result of a multiplication between two ('2ModPMap'*'2ModPMap' -> '2ModPMult')
%       algType - (String) Name of the used random walk algorithm <'probRand' | 'probRandExt' | ... >
%       opLog   - Structure containing the operation history as structure
%
%
%   return values
%     res:  the probStruct
%     errStr: If an error occured, errStr is identify the error
%
%
% Bjoern W. Kreher
% 5/08
%
% UNIX


res= [];   errStr= '';

if nargin < 1
    errStr= 'probstruct_init(error): Arg. is not of type mr- or probStruct';
    return;   
end

% seems already a probmap, which have to be converted
[ok, errStr, verStr]= probstruct_istype(mrStruct);
if (nargin == 1)
    if strcmp(verStr, 'V1')
        [res, errStr]= local_convertV1_to_V2(mrStruct);
        return;
    elseif strcmp(verStr, 'V2')
        res= mrStruct;
        return
    elseif mrstruct_istype(mrStruct)
        [res, errStr]= local_genDummy(mrStruct);
        return;        
    end
end

if nargin < 5
    errStr= 'probstruct_init(error): not enough arguments';
    return;   
end
   
if ~mrstruct_istype(mrStruct)
    errStr= 'probstruct_init(error): Arg. is not of type mr- or probStruct';
    return;   
end


[pInfo, errStr]= local_getBaseStruct(mapName, dataType, algType, opLog);
mrStruct.user.probInfo= pInfo;
res= mrStruct;





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
function [pInfo, errStr]= local_getBaseStruct(mapName, dataType, algType, opLog)
pInfo= []; errStr= '';

dTCell= {'2ModMult'; 'Mult'; '2ModPMap'; 'PMap'};
if ~any(strcmp(dTCell, dataType))
    errStr= sprintf('%s(local_getBaseStruct): ''%s'' is a non valid dataType', ...
        mfilename, dataType);
    return;
end

    pInfo.mapName       = mapName;      %'li_V1*li_CGL'
    pInfo.algType       = algType;      %'ProbRandExt'
    pInfo.dataType      = dataType;     %'2ModMult'
    pInfo.operationLog  = opLog;
%    pInfo.operationLog.type= '';
%    pInfo.operationLog.date= '';
%    pInfo.operationLog.params= opLog;
    pInfo.histStruc     = [];
    pInfo.logHistStruc  = [];
    pInfo.MinMaxStruc   = [];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getBaseStruct
%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [res, errStr]= local_genDummy(mrStruct)
res= []; errStr= '';

sizeAy= mrstruct_query(mrStruct, 'sizeAy');

if length(sizeAy) ~= 3
    errStr= sprintf('%s::local_genDummy(error): mrStruct has wrong dimension', mfilename);
    return;
end

%% generate opLog
opLog.type= 'import';
opLog.date= date;
opLog.params= {};

[res, errStr]= probstruct_init(mrStruct, 'import', 'PMap', 'undef', opLog);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [res, errStr]= local_convertV1_to_V2(pStructV1)
res= []; errStr= '';

%%  save and clear old user entry
oStrc= pStructV1.user;
pStructV1.user= rmfield(pStructV1.user, 'mapName');
pStructV1.user= rmfield(pStructV1.user, 'algName');
pStructV1.user= rmfield(pStructV1.user, 'versionLog');
pStructV1.user= rmfield(pStructV1.user, 'genDate');
pStructV1.user= rmfield(pStructV1.user, 'operationLog');

%% extract probAlgName
idx= strfind(oStrc.algName, '_V');
if isempty(idx)
    errStr= sprintf('%s(local_convertV1_to_V2): ''%s'' is a non valid algorith name', ...
        mfilename, oStrc.algName);
    return;
end
algType= oStrc.algName(1:(idx(end) - 1));

%% check if already multiplied
if ~isempty(strfind(oStrc.operationLog.operationLog, ' * '))
    if strcmp(algType, 'ProbRandExt')
        dataType= '2ModMult';        
    else
        dataType= 'Mult';        
    end
else
    if strcmp(algType, 'ProbRandExt')
        dataType= '2ModPMap';        
    else
        dataType= 'PMap';        
    end
end

%% generate opLog
opLog.type= 'undef';
opLog.date= oStrc.genDate;
opLog.params= oStrc.operationLog;
opLog.params.algName= oStrc.algName;


[res, errStr]= probstruct_init(pStructV1, oStrc.mapName, dataType, algType, opLog);
