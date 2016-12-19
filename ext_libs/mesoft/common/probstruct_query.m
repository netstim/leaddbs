function [res, errStr, res2]= probstruct_query(varargin)
%extracts information from a probStruct concerning the probMaps. Futher information (spatial) can be achived by mrstruct_query 
% function [res, errStr, res2]= probstruct_query(probStruct, commandStr[, param1[, param2[... [, paramN]]]])
%
%   probStruct: should contain a valid probStruct
%   command:    {'probType' | 'mapTypes' | 'mapName' | 'getMap' | 'getMinMax' | 'getHist' | 'getLogHist' | 'cmpDimensions'}
%   param[1..end] the parameter specified by the command string
% 
% 'probType'     - returns the probType of the probStruct:
%                  'PMap' Standard visiting map. Not multipied and one value for each voxel
%                  'Mult' Product of two or mor 'PMaps'
%                  '2ModPMap' Visiting maps conatining also directional
%                     information (two values for each voxel) (cmp. probMap paper NeuroImage)
%                  '2ModMult' Specialized product of two '2ModPMap' (cmp. probMap paper NeuroImage)
% 'mapTypes'     - returns a list of names of map types the current probStruct does support:
%                  'PMap': supports 'probMap'
%                  'Mult': supports 'probMap'
%                  '2ModPMap': supports 'probMap' corresponds to the sum of both maps
%                  '2ModMult': supports 'sumMap' Probability Index of a connection to both seed regions (PIBS); 
%                              supports 'conMap' Probability Index of forming part of BOI (PIBI); 
%                              supports 'mergMap' PIBS - PIBI; 
%                              supports 'conWeight' Probability Index of connecting fibre configuration (PICC)
%                              supports 'conFrac' Fraction of connecting fibre configuration (FCC) is not supportd by all maps
%                              (cmp. probMap paper NeuroImage)
% 'mapName'      - returns the name of the probStruct
% 'getMap'       - returns a data array containing the data of mapType defined in param1 (string)
% 'getMinMax'    - returns the min (res(1)) and max (res(2)) value of the data array containing the values of the 
%                  mapType defined in param1 (string). 
%                  If the probstruct_modify 'createHists' command was applied, this procedure is much faster
% 'getHist'      - returns the normalized histogram of the the data of mapType defined in param1 (string). 
%                  res(1, :) contains the bins position, res(2, :) contains the fraction.
%                  If the probstruct_modify 'createHists' command was applied, this procedure is much faster
% 'getLogHist'   - returns the normalized log10 histogram of the the data of mapType defined in param1 (string). 
%                  res(1, :) contains the bins position, res(2, :) contains the fraction.
%                  If the probstruct_modify 'createHists' command was applied, this procedure is much faster
% 'cmpDimensions'- Compares, if data fits to the spacial properties of a given dtd/mrStruct (param1)
%                  returns < 0: data doesn't fit and it is not prossible to reslice
%                  returns = 0: data doesn't fit but reslicing is possible
%                  returns > 1: data fit and no reslicing is necessary
%
% return values:
%   res: The result of the method or if an error occured res is empty
%   errStr: a message string, specifying the error which occured
%   res2: additional output argument, depending on the command
% 
% Description of the different commands: 
%
% Bjoern W. Kreher
% 01/08
%
% UNIX


res= [];    res2= [];    errStr= '';
%% check input param
if length(varargin) < 2
    errStr= sprintf('%s(error): There have to be at least two parameters', mfilename);
    return;
end

if (~isempty(varargin)) && probstruct_istype(varargin{1})
    probStruct= varargin{1};
else
    errStr= sprintf('%s(error): First param have to be probStruct', mfilename);
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
if strcmp(commandStr, 'probType')
    [res, errStr]= local_getProbType(probStruct);
elseif strcmp(commandStr, 'mapTypes')
    [res, errStr]= local_getMapTypes(probStruct);
elseif strcmp(commandStr, 'mapName')
    [res, errStr]= local_getMapName(probStruct);
elseif strcmp(commandStr, 'getMap')
    if ischar(param{1})
        mapTypeStr= param{1};
    else
        errStr= sprintf('%s(getMap): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_getMap(probStruct, mapTypeStr);
elseif strcmp(commandStr, 'getMinMax')
    if ischar(param{1})
        mapTypeStr= param{1};
    else
        errStr= sprintf('%s(getMinMax): invalid data', mfilename);
        return
    end    
    [res, errStr]= local_getMinMax(probStruct, mapTypeStr);
elseif strcmp(commandStr, 'getHist')
    if ischar(param{1})
        mapTypeStr= param{1};
    else
        errStr= sprintf('%s(getMap): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_getHist(probStruct, mapTypeStr, 'lin', 100);
elseif strcmp(commandStr, 'getLogHist')
    if ischar(param{1})
        mapTypeStr= param{1};
    else
        errStr= sprintf('%s(getMap): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_getHist(probStruct, mapTypeStr, 'log', 100);
elseif strcmp(commandStr, 'cmpDimensions')
    if mrstruct_istype(param{1})
        mrStruct= param{1};
    else
        errStr= sprintf('%s(cmpDimensions): invalid data', mfilename);
        return
    end    
    [res, errStr]= mrstruct_query(probStruct, 'cmpDimensions', mrStruct);    
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
function [res, errStr]= local_getMapName(pStruct)
res= []; errStr= '';

res= pStruct.user.probInfo.mapName;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getProbType(pStruct)
res= []; errStr= '';

res= pStruct.user.probInfo.dataType;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getMapTypes(pStruct)
errStr= '';
res= [];


probType= local_getProbType(pStruct);

if strcmp(probType, 'PMap') || strcmp(probType, 'Mult')
    res= {'probMap'};
elseif strcmp(probType, '2ModMult')
    res= {'sumMap'; 'conMap'; 'mergMap'; 'conWeight'};
    
    if strcmp(pStruct.user.probInfo.operationLog.type, 'conMULT') && ...
            (length(pStruct.user.probInfo.operationLog.params) == 3) && ...
            isfield(pStruct.user.probInfo.operationLog.params{3}, 'conFractionAy')
        res(end+1)= {'conFrac'};
    end            
elseif strcmp(probType, '2ModPMap')
    res= {'sumMap'};
else
    errStr= sprintf('%s(local_getMapTypes): ''%s'' is not a supported probMap type', ...
            mfilename, probType);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, newType]= local_getMap(pStruct, mapType)
res= []; errStr= ''; newType= '';


if ~any(strcmp(local_getMapTypes(pStruct), mapType))
    errStr= sprintf('%s(local_getMap): ''%s'' is not a supported by the given probStruct', ...
        mfilename, mapType);
    return;
end

if strcmp(mapType, 'probMap')
    res= pStruct.dataAy; newType= pStruct.user.probInfo.dataType;
elseif strcmp(mapType, 'sumMap')
    if strcmp(pStruct.user.probInfo.dataType, '2ModPMap')
        newType= 'PMap';
    elseif strcmp(pStruct.user.probInfo.dataType, '2ModMult')
        newType= 'Mult';
    else
        errStr= sprintf('%s(local_getMap): dataType ''%s'' is not a supported yet', ...
            mfilename, pStruct.user.probInfo.dataType);
        return;
    end
    res= sum(pStruct.dataAy, 4);    
elseif strcmp(mapType, 'conMap')
    newType= 'Mult';
    res= pStruct.dataAy(:, :, :, 1);
elseif strcmp(mapType, 'mergMap')
    newType= 'Mult';
    res= pStruct.dataAy(:, :, :, 2);
elseif strcmp(mapType, 'conWeight')
    newType= 'Mult';
    sizeAy= mrstruct_query(pStruct, 'sizeAy');
    res= NaN(sizeAy(1:3));
    ss= sum(pStruct.dataAy, 4);
    idx= find(ss > 0);
    res(idx)= pStruct.dataAy(idx)./ss(idx);
elseif strcmp(mapType, 'conFrac')
        if strcmp(pStruct.user.probInfo.operationLog.type, 'conMULT') && ...
            (length(pStruct.user.probInfo.operationLog.params) == 3) && ...
            isfield(pStruct.user.probInfo.operationLog.params{3}, 'conFractionAy')
            res= pStruct.user.probInfo.operationLog.params{3}.conFractionAy;
        else
            errStr= sprintf('%s(local_getMap): mapType ''conFrac'' is not a supported for this kind of probStruct', mfilename);
            return;
        end
else
    errStr= sprintf('%s(local_getMap): ''%s'' should, but is not supported', ...
        mfilename, mapType);
    return;    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getMinMax(pStruct, mapType)
res= []; errStr= '';

if isfield(pStruct.user.probInfo.MinMaxStruc, mapType)
        res= pStruct.user.probInfo.MinMaxStruc.(mapType);
        return
end

[map, errStr]= local_getMap(pStruct, mapType);
if isempty(map)
    return
end
res= [min(map(:)) max(map(:))];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, pAy]= local_getHist(pStruct, mapType, linOrLog, binNo)
res= []; errStr= ''; pAy= [];

if strcmp(linOrLog, 'lin')
    if isfield(pStruct.user.probInfo.histStruc, mapType)
        res= pStruct.user.probInfo.histStruc.(mapType);
        return
    end
elseif strcmp(linOrLog, 'log')
    if isfield(pStruct.user.probInfo.logHistStruc, mapType)
        res= pStruct.user.probInfo.logHistStruc.(mapType);
        return
    end
else
    errStr= sprintf('%s(local_getMap): scale type ''%s'' is not supported yet', ...
        mfilename, linOrLog);
    return;        
end

[map, errStr]= local_getMap(pStruct, mapType);
if isempty(map)
    return
end

idx= find(map > 0);
if isempty(idx)
    res= [0 1; 0 0];
    pAy= [0 0.5 1];
    errStr= sprintf('%s(local_getMap): map is empty', mfilename);
    return
end
if strcmp(linOrLog, 'log')
    vals= log10(map(idx));
    pAy= 10.^linspace(min(vals), max(vals), binNo + 1);
    [A, B]= hist(log10(map(idx)), binNo);
    res= [10.^B; A./diff(pAy)];
else
    vals= map(idx);
    pAy= linspace(min(vals), max(vals), binNo + 1);
    [A, B]= hist(map(idx), binNo);
    res= [B; A./diff(pAy)];
end



