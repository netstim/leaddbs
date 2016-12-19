function [res, errStr, res2]= probstruct_modify(varargin)
%methods to change properties of a probstruct
%   function [res, errStr, res2]= probstruct_modify(probStruct, command, parm1[, parm2[, ...[, parmN]]])
%
%   probStruct: A valid probstruct, which will be modified
%   command:   spezifies the operation which should be applied
%     {'createHists' | 'mapName'}
%   param[1..end] the parameter specified by the command string
%
%   'createHists' - Creates for histograms and min and max values for a faster access using 
%                   the probstruct_query command 'getHist', 'getLogHist', or 'getMinMax'.
%                   The createHist command calcualtes the values and saves it inside the probStruct 
%                   returned by this command.
%   'mapName'     - Set the mapName of the probStruct to the string param1
% 
% return values:
%   return: the modified probStruct or if an error occured the old one
%   errStr: a message string, specifying the error which occured
%   res2: additional output argument, depending on the command
% 
%
%  B.W. Kreher
%  01/08
%
%  UNIX
%

res= [];    res2= [];   errStr= '';
%% check input param
if nargin < 2
    errStr= sprintf('%s(error): There have to be at least two parameters', mfilename);
    return;
end

[ok, errStr, verStr]= probstruct_istype(varargin{1});
if ok 
    if strcmp(verStr, 'V2')
        pStruct= varargin{1};
    else
        [pStruct, errStr]= probstruct_init(varargin{1});
        if isempty(pStruct)
            return
        end
    end
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

%% transfer input param
paramMax= 10;
param= cell(paramMax, 1);
for i= 1:paramMax
    if length(varargin) < (i + 2)
        param{i}= [];
    else
        param{i}= varargin{i + 2};
    end
end

res= pStruct;

%% split to different commands
if strcmp(commandStr, 'createHists')
    [res, errStr]= local_createHists(pStruct);
elseif strcmp(commandStr, 'mapName')
    if ischar(param{1})
        mapName= param{1};
    else
        errStr= sprintf('%s(mapName): invalid data', mfilename);
        return
    end    

    [res, errStr]= local_setMapName(pStruct, mapName);
else
    errStr= sprintf('%s(error): Command ''%s'' is not implemented', mfilename, commandStr);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pStruct, errStr]= local_setMapName(pStruct, mapName)
res= []; errStr= '';

pStruct.user.probInfo.mapName= mapName;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pStruct, errStr]= local_createHists(pStruct)

errStr= '';
pStruct;

[mapStrs, errStr]= probstruct_query(pStruct, 'mapTypes');
if isempty(mapStrs)
    return
end

histStruc = []; 
logHistStruc= []; 
MinMaxStruc= [];

for i= 1:length(mapStrs)
   [histAy, errStr, pAy]= probstruct_query(pStruct, 'getLogHist', mapStrs{i});
   if isempty(histAy)
        return
   end    
   logHistStruc.(mapStrs{i}) = histAy;
   
   [histAy, errStr, pAy]= probstruct_query(pStruct, 'getHist', mapStrs{i});
   if isempty(histAy)
        return
   end    
   histStruc.(mapStrs{i}) = histAy;
   MinMaxStruc.(mapStrs{i}) = [min(pAy) max(pAy)];
end   
   
pStruct.user.probInfo.histStruc= histStruc;
pStruct.user.probInfo.logHistStruc= logHistStruc;
pStruct.user.probInfo.MinMaxStruc= MinMaxStruc;
