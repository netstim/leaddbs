function [res, errStr, res2]= probstruct_op(varargin)
%special operation between and on probStuct, including the random walk methods
%   function [res, errStr, res2]= probstruct_op(command, parm1[, parm2[, ...[, parmN]]])
%
%   command:   spezifies the operation which should be applied
%     {'ADD' | 'MULT' | 'NORM' | 'fProbAlgs' | 'fProbAlgArgs' | 'extract_???' | 'apply_???'}
%   param[1..end] the parameter specified by the command string
%
% 'ADD'         - Add several probStruct (param1, param2, ...)with the same mapType together
% 'MULT'        - Multiplies several (param1, param2, ...) probStruct pointwise.
%                 rules: PMap * PMap -> Mult; Mult * Mult -> Mult
%                       2ModPMap X 2ModPMap -> 2ModMult (The X remarks the connection weighted multiplication)
%                       2ModMult * 2ModMult -> 2ModMult
% 'NORM'        - Normalizes a probStruct (param1) to get a maximum of one. In case of 2ModXXX maps,
%                 the max of the sum of both maps is used to determine the  
%                 normalization factor
% 'extract_???' - extracts a map from the probStruct (param1) and returns it as probStruct. 
%                 The possible map names(???) which is supported by the probStruct results from the
%                       probstruct_query(probStruct, 'mapTypes')
% 'fProbAlgs'   - search the library for available random walk methods and
%                 returns a list
% 'fProbAlgArgs'- For the specified randomwalk method (param1), the method
%                 returns a list of the paramter names (res) and a list of corrsponding
%                 default values for the parameter. 
% 
% 'apply_???'   - applies the specified RW algorithm (???, see fProbAlgs). 
%                 param1: dtdStruct
%                 param2: seedVc {Nx3 int} specifying the different seed voxel
%                 param3: maskAy {size(dtsStruct) bool} mask specifying the
%                         trackable are (eg. WM mask) 
%                 param4: cellarray containing the alg specific parameters
%                         (see also fProbAlgArgs)
%                 param5: handle for a GUI element for verbose messages
% 
%
% return values:
%   return: a probStruct or if an error occured the old one or just the
%           empyt set
%   errStr: a message string, specifying the error which occured
%   res2: additional output argument, depending on the command
% 
%
%  Bjoern W. Kreher
%  01/08
%
%  UNIX
%

res= [];    res2= [];   errStr= '';
%% check input param
if nargin < 1
    errStr= sprintf('%s(error): There have to be at least one parameters', mfilename);
    return;
end

if ~ischar(varargin{1})
    errStr= sprintf('%s(error): Command have to be a string', mfilename);
    return;
else
    commandStr= varargin{1};
end

%% transfer input param
probArgNo= 0;
argNo= 0;
paramMax= 10;
param= cell(paramMax, 1);
for i= 1:paramMax
    if length(varargin) < (i + 1)
        param{i}= [];
        if probArgNo == 0   % save the number of probmaps 
            probArgNo= i - 1;
        end
    else
        if ~isempty(varargin{i + 1})
            argNo= i;
        end
        if probstruct_istype(varargin{i + 1}) % be sure, that probstruct is of V2
            param{i}= probstruct_init(varargin{i + 1});
        else
            param{i}= varargin{i + 1};
            if probArgNo == 0   % save the number of probmaps 
                probArgNo= i - 1;
            end
        end
    end
end

extractStr= 'extract_';
randomWalkStr= 'apply_';
%% split to different commands
if strcmp(commandStr, 'fProbAlgs')
    [res, errStr, res2]= local_fProbAlgs;
elseif strcmp(commandStr, 'fProbAlgArgs')
    if ~isempty(param{1}) && ischar(param{1})
        algName= param{1};
    else
        errStr= sprintf('%s(fProbAlgParam): Argument have to be a string', mfilename);
        return
    end
    [res, errStr, res2]= local_fProbAlgArgs(algName);
elseif strcmp(commandStr, 'ADD')
    if (probArgNo >= 2) && (isempty(param{probArgNo + 1}) || ischar(param{probArgNo + 1}))
        pStructCell= param(1:probArgNo); mapName= param{probArgNo + 1};
    else
        errStr= sprintf('%s(ADD): Two or more probStructs have to be defined', mfilename);
        return
    end
    [res, errStr]= local_add(pStructCell, mapName);
elseif strcmp(commandStr, 'MULT')
    if (probArgNo >= 2) && (isempty(param{probArgNo + 1}) || ischar(param{probArgNo + 1})) && (isempty(param{probArgNo + 2}) || isscalar(param{probArgNo + 2})) 
        pStructCell= param(1:probArgNo); mapName= param{probArgNo + 1}; B= param{probArgNo + 2};
    else
        errStr= sprintf('%s(MULT): Two or more probStructs have to be defined', mfilename);
        return
    end
    [res, errStr]= local_mult(pStructCell, mapName, B);
elseif strcmp(commandStr, 'NORM')
    if probstruct_istype(param{1}) && (isempty(param{2}) || ischar(param{2}))
        pStruct= param{1}; mapName= param{2};
    else
        errStr= sprintf('%s(NORM): Not expected parameter', mfilename);
        return
    end
    [res, errStr]= local_norm(pStruct, mapName);
elseif isequal(strfind(commandStr, extractStr), 1)
    if probstruct_istype(param{1}) && (isempty(param{2}) || ischar(param{2}))
        pStruct= param{1}; mapName= param{2}; 
        mapType= commandStr((length(extractStr)+1):end);
    else
        errStr= sprintf('%s(extract_XXX): Not expected parameter', mfilename);
        return
    end
    [res, errStr]= local_extract(pStruct, mapName, mapType);

elseif isequal(strfind(commandStr, randomWalkStr), 1)
    if dtdstruct_istype(param{1}) && ...
            isnumeric(param{2}) && (size(param{2}, 2) == 3) && isequal(param{2}, round(param{2})) && ...
            islogical(param{3})
        % last argument is handle for status line (if handle)
        if ishandle(param{argNo})
            vHd= param{argNo}; 
            argNo= argNo - 1;
        else
            vHd= [];
        end
        
        algName= commandStr((length(randomWalkStr)+1):end);
        dtdStruct= param{1}; seedVc= param{2}; maskAy= param{3}; 
        trackParams= param(4:argNo);
    else
        errStr= sprintf('%s(%sXXX): Not expected parameter', mfilename, randomWalkStr);
        return
    end
    [res, errStr]= local_startRWAlg(algName, dtdStruct, seedVc, maskAy, trackParams, vHd);

else
    errStr= sprintf('%s(error): Command ''%s'' is not implemented', mfilename, commandStr);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, algVers]= local_fProbAlgs(algName)
res= {}; errStr= ''; algVers= {};

algNameCell      = {'ProbRand'; 'ProbRandExt'; 'ProbRandExtHARDI';'ProbRandRef'; 'ProbFACT'; 'ProbKoch'; 'ProbTarget'};
dllNameCell      = {'random_walk_prob'; 'random_walk_prob';'random_walk_prob'; 'random_walk_ref'; 'prob_fakt_main'; 'koch_alg_c'; 'random_walk_goal'};

for i= 1:length(dllNameCell)
    if exist(dllNameCell{i}, 'file') == 3
        res{end + 1, 1}= algNameCell{i};
        verStr= private_catchAlgVersionName(algNameCell{i});
        algVers{end + 1, 1}= verStr;
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, defaultParam]= local_fProbAlgArgs(algName)
res= []; errStr= ''; defaultParam= [];

algNameCell      = {'ProbRand'; 'ProbRandExt'; 'ProbRandExtHARDI'; 'ProbRandRef'; 'ProbFACT'; 'ProbKoch'; 'ProbTarget'};
algParamCell     = {'FiberLen:', 'Walk No:', 'exp A:', '', ''; ...
                    'FiberLen:', 'Walk No:', 'exp A:', 'Allow Revisits', ''; ...
                    'FiberLen:', 'Walk No:', 'exp A:', 'Allow Revisits', 'HARDIdata' ;...
                    'FiberLen:', 'Walk No:', 'exp A:', '', ''; ...
                    'FiberLen:', 'Sample P:', '', '', ''; ...
                    'JumpNo:', 'Walk No:', 'exp A:', '', ''; ...
                    'FiberLen:', 'Walk No:', 'exp A:', '', '';};
param1           = [ 150   150 150 5000 1000  100  150];
param2           = [1000  1000 1000 200   10 1000 1000];
param3           = [   2     4  4   2  NaN    7    4];
param4           = [ NaN     1  1  NaN  NaN  NaN  NaN];

idx= find(strcmp(algNameCell, algName));
res= algParamCell(idx, :);
defaultParam= [param1(idx), param2(idx), param3(idx), param4(idx)];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_norm(pStruct, mapName)
res= []; errStr= '';

reS= pStruct;
% reS.user.probInfo.algType= ; %bleibt
% reS.user.probInfo.dataType= ; %bleibt
reS.user.probInfo.logHistStruc= [];
reS.user.probInfo.histStruc= [];
reS.user.probInfo.MinMaxStruc= [];
res.user.probInfo.operationLog= [];
reS.user.probInfo.operationLog.type= 'NORM';
reS.user.probInfo.operationLog.date= datestr(now);
reS.user.probInfo.operationLog.params{1}= pStruct.user.probInfo.operationLog;

% generate new mapName
mName= sprintf('NORM(%s)', pStruct.user.probInfo.mapName);

if any(strcmp(pStruct.user.probInfo.dataType, {'2ModMult'; '2ModPMap'}))
    tmpAy= sum(pStruct.dataAy, 4);
    reS.dataAy= pStruct.dataAy/max(tmpAy(:));
elseif any(strcmp(pStruct.user.probInfo.dataType, {'Mult'; 'PMap'}))
    reS.dataAy= pStruct.dataAy/max(pStruct.dataAy(:));
else
    errStr= sprintf('%s::local_norm(error): The dataRype ''%s'' is not supported yet', ...
        mfilename, pStruct.user.probInfo.dataType);
    return
end
mapName= mName;
reS.user.probInfo.mapName= mapName;
res= reS;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_add(pStructCell, mapName)
res= []; errStr= '';

% init result struct
reS= pStructCell{1};
% reS.user.probInfo.algType= ; %bleibt
% reS.user.probInfo.dataType= ; %bleibt
reS.user.probInfo.logHistStruc= [];
reS.user.probInfo.histStruc= [];
reS.user.probInfo.MinMaxStruc= [];
reS.user.probInfo.operationLog= [];
reS.user.probInfo.operationLog.type= 'ADD';
reS.user.probInfo.operationLog.date= datestr(now);
reS.user.probInfo.operationLog.params= cell(length(pStructCell), 1);
reS.user.probInfo.operationLog.params{1}= pStructCell{1}.user.probInfo.operationLog;

% generate new mapName
mName= sprintf('(%s)', pStructCell{1}.user.probInfo.mapName);

for i= 2:length(pStructCell)
    % TEst data compatibility
    if ~strcmp(reS.user.probInfo.algType, pStructCell{i}.user.probInfo.algType)
        errStr= sprintf('%s::local_add(error): The maps origin from different methods', mfilename);
        return
    end
    if ~strcmp(reS.user.probInfo.dataType, pStructCell{i}.user.probInfo.dataType)
        errStr= sprintf('%s::local_add(error): The maps have different types', mfilename);
        return
    end
    [ok, errStr]= probstruct_query(reS, 'cmpDimensions', pStructCell{i});
    if (ok ~= 2), return; end
    
    % Add map together
    reS.dataAy= reS.dataAy + pStructCell{i}.dataAy;
    reS.user.probInfo.operationLog.params{i}= pStructCell{i}.user.probInfo.operationLog;

    % generate new mapName
    mName= sprintf('%s + (%s)', mName, pStructCell{i}.user.probInfo.mapName);
end

if isempty(mapName)
    mapName= mName;
end
reS.user.probInfo.mapName= mapName;
res= reS;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_mult(pStructCell, mapName, B)
res= []; errStr= '';

if isempty(B)
    B= 0.05;
end
% init result struct
reS= pStructCell{1};
% reS.user.probInfo.algType= ; %bleibt
% reS.user.probInfo.dataType= ; %bleibt
reS.user.probInfo.logHistStruc= [];
reS.user.probInfo.histStruc= [];
reS.user.probInfo.MinMaxStruc= [];
reS.user.probInfo.operationLog= [];
reS.user.probInfo.operationLog.type= 'MULT';
reS.user.probInfo.operationLog.date= datestr(now);
reS.user.probInfo.operationLog.params= cell(length(pStructCell), 1);
reS.user.probInfo.operationLog.params{1}= pStructCell{1}.user.probInfo.operationLog;

if strcmp(pStructCell{1}.user.probInfo.dataType, '2ModPMap')
    reS.user.probInfo.dataType= '2ModMult';
elseif strcmp(pStructCell{1}.user.probInfo.dataType, 'PMap')
    reS.user.probInfo.dataType= 'Mult';
else
    reS.user.probInfo.dataType= pStructCell{1}.user.probInfo.dataType;
end
    
% generate new mapName
mName= sprintf('(%s)', pStructCell{1}.user.probInfo.mapName);


for i= 2:length(pStructCell)
    % TEst data compatibility
    if ~strcmp(pStructCell{1}.user.probInfo.algType, pStructCell{i}.user.probInfo.algType)
        errStr= sprintf('%s::local_mult(error): The maps origin from different methods', mfilename);
        return
    end
    if ~strcmp(pStructCell{1}.user.probInfo.dataType, pStructCell{i}.user.probInfo.dataType)
        errStr= sprintf('%s::local_mult(error): The maps have different types', mfilename);
        return
    end
    [ok, errStr]= probstruct_query(reS, 'cmpDimensions', pStructCell{i});
    if (ok ~= 2), return; end
    
    % Multiply map together
    reS.user.probInfo.operationLog.params{i}= pStructCell{i}.user.probInfo.operationLog;
    if strcmp(pStructCell{i}.user.probInfo.dataType, '2ModPMap') % connected weighted multiplication
        if i > 2
            errStr= sprintf('%s::local_mult(error): For dataType ''2ModPMap'' only two candidates are allowed', mfilename);
            return
        end

        [reS.dataAy, errStr, fracAy]= private_connWeightMult(pStructCell{1}.dataAy, pStructCell{2}.dataAy, B);
        reS.user.probInfo.operationLog.type= 'conMULT';
        reS.user.probInfo.operationLog.params{i+1}.exp_B= B;
        reS.user.probInfo.operationLog.params{i+1}.conFractionAy= fracAy;
        
        mName= sprintf('%s X (%s)', mName, pStructCell{i}.user.probInfo.mapName);        
    else % pointwise multiplication
        reS.dataAy= reS.dataAy .* pStructCell{i}.dataAy;
        mName= sprintf('%s * (%s)', mName, pStructCell{i}.user.probInfo.mapName);
    end

    % generate new mapName
end

if isempty(mapName)
    mapName= mName;
end
reS.user.probInfo.mapName= mapName;
res= reS;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_extract(pStruct, mName, mapType)
errStr= '';
res= [];

res= mrstruct_init('volume', [], pStruct);
[res.dataAy, errStr, dataType]= probstruct_query(pStruct, 'getMap', mapType);
if isempty(res.dataAy),  res= []; return; end

res.user.probInfo.dataType= dataType;
%res.user.probInfo.algType= ; %bleibt
res.user.probInfo.logHistStruc= [];
res.user.probInfo.histStruc= [];
res.user.probInfo.MinMaxStruc= [];
res.user.probInfo.operationLog= [];
res.user.probInfo.operationLog.type= mapType;
res.user.probInfo.operationLog.date= datestr(now);
res.user.probInfo.operationLog.params{1}= pStruct.user.probInfo.operationLog;

% generate new mapName
if isempty(mName)
    res.user.probInfo.mapName= sprintf('%s(%s)', mapType, pStruct.user.probInfo.mapName);
else
    res.user.probInfo.mapName= mName;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_startRWAlg(algName, dtdStruct, seedVc, maskAy, trackParams, vHd)
res= []; errStr= ''; 

% Hole tensor daten
eigVec= dtdstruct_query(dtdStruct, 'getEigVec');
eigVal= dtdstruct_query(dtdStruct, 'getEigVal');
sizeAy= dtdstruct_query(dtdStruct, 'sizeAy');
vox= dtdstruct_query(dtdStruct, 'vox');

% teste allgemeine parameter
if isempty(local_fProbAlgArgs(algName))
    errStr= sprintf('%s::local_startRWAlg(error): algorithm is not supported here', mfilename);
    return
end 
if ~isequal(sizeAy, size(maskAy))
    errStr= sprintf('%s::local_startRWAlg(error): maskAy has wrong size', mfilename);
    return
end 
for i= 1:3
    if (min(seedVc(:, i)) < 1) || (max(seedVc(:, i)) > sizeAy(i))
        errStr= sprintf('%s::local_startRWAlg(error): seed points are out of volume', mfilename);
        return
    end         
end

% Bereite spezielle parameter vor und starte RW
if strcmp(algName, 'ProbRand')
    if length(trackParams) < 3
        errStr= sprintf('%s::local_startRWAlg::ProbRandExt(error): not enough parameters', mfilename);
        return
    end
    if isscalar(trackParams{1}) && (trackParams{1} == abs(round(trackParams{1})))
        fibMaxLen= trackParams{1};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbRand(error): fibMaxLen has wrong properties', mfilename);
        return
    end
    if isscalar(trackParams{2}) && (trackParams{2} == abs(round(trackParams{2})))
        repNo= trackParams{2};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbRand(error): repNo has wrong properties', mfilename);
        return
    end
    if isscalar(trackParams{3}) && (trackParams{3} == abs(trackParams{3}))
        exp_alpha= trackParams{3};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbRand(error): exp_alpha has wrong properties', mfilename);
        return
    end
    [mapAy, errStr, opLog, dataType]= private_ProbRand(eigVec.dataAy, eigVal.dataAy, vox, seedVc, maskAy, ...
        repNo, fibMaxLen, exp_alpha, vHd);
    if isempty(mapAy), return; end
elseif strcmp(algName, 'ProbRandExt')
    if length(trackParams) < 4
        errStr= sprintf('%s::local_startRWAlg::ProbRandExt(error): not enough parameters', mfilename);
        return
    end
    if isscalar(trackParams{1}) && (trackParams{1} == abs(round(trackParams{1})))
        fibMaxLen= trackParams{1};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbRandExt(error): fibMaxLen has wrong properties', mfilename);
        return
    end
    if isscalar(trackParams{2}) && (trackParams{2} == abs(round(trackParams{2})))
        repNo= trackParams{2};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbRandExt(error): repNo has wrong properties', mfilename);
        return
    end
    if isscalar(trackParams{3}) && (trackParams{3} == abs(trackParams{3}))
        exp_alpha= trackParams{3};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbRandExt(error): exp_alpha has wrong properties', mfilename);
        return
    end
    if isscalar(trackParams{4}) 
        allowRevisit= trackParams{4} == 1;
    else
        errStr= sprintf('%s::local_startRWAlg::ProbRandExt(error): exp_alpha has wrong properties', mfilename);
        return
    end
    [mapAy, errStr, opLog, dataType]= private_ProbRandExt(eigVec.dataAy, eigVal.dataAy, vox, seedVc, maskAy, ...
        repNo, fibMaxLen, exp_alpha, allowRevisit, vHd);
    if isempty(mapAy), return; end
elseif strcmp(algName, 'ProbRandExtHARDI')
    if length(trackParams) < 4
        errStr= sprintf('%s::local_startRWAlg::ProbRandExt(error): not enough parameters', mfilename);
        return
    end
    if isscalar(trackParams{1}) && (trackParams{1} == abs(round(trackParams{1})))
        fibMaxLen= trackParams{1};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbRandExt(error): fibMaxLen has wrong properties', mfilename);
        return
    end
    if isscalar(trackParams{2}) && (trackParams{2} == abs(round(trackParams{2})))
        repNo= trackParams{2};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbRandExt(error): repNo has wrong properties', mfilename);
        return
    end
    if isscalar(trackParams{3}) && (trackParams{3} == abs(trackParams{3}))
        exp_alpha= trackParams{3};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbRandExt(error): exp_alpha has wrong properties', mfilename);
        return
    end
    if isscalar(trackParams{4}) 
        allowRevisit= trackParams{4} == 1;
    else
        errStr= sprintf('%s::local_startRWAlg::ProbRandExt(error): exp_alpha has wrong properties', mfilename);
        return
    end
    if true,
        hardiData = trackParams{5},
    else
        errStr= sprintf('%s::local_startRWAlg::ProbRandExt(error): hardi data not compatible', mfilename);
        return
    end
    
    [mapAy, errStr, opLog, dataType]= private_ProbRandExtHARDI(eigVec.dataAy, eigVal.dataAy,hardiData, vox, seedVc, maskAy, ...
        repNo, fibMaxLen, exp_alpha, allowRevisit, vHd);
    if isempty(mapAy), return; end    
elseif strcmp(algName, 'ProbRandRef')
    if length(trackParams) < 3
        errStr= sprintf('%s::local_startRWAlg::ProbRandRef(error): not enough parameters', mfilename);
        return
    end
    if isscalar(trackParams{1}) && (trackParams{1} == abs(round(trackParams{1})))
        fibMaxLen= trackParams{1};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbRandRef(error): fibMaxLen has wrong properties', mfilename);
        return
    end
    if isscalar(trackParams{2}) && (trackParams{2} == abs(round(trackParams{2})))
        repNo= trackParams{2};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbRandRef(error): repNo has wrong properties', mfilename);
        return
    end
    if isscalar(trackParams{3}) && (trackParams{3} == abs(trackParams{3}))
        exp_alpha= trackParams{3};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbRandRef(error): exp_alpha has wrong properties', mfilename);
        return
    end
    [mapAy, errStr, opLog, dataType]= private_ProbRandRef(eigVec.dataAy, eigVal.dataAy, vox, seedVc, maskAy, ...
        repNo, fibMaxLen, exp_alpha, vHd);
    if isempty(mapAy), return; end
elseif strcmp(algName, 'ProbFACT')
    if length(trackParams) < 3
        errStr= sprintf('%s::local_startRWAlg::ProbFACT(error): not enough parameters', mfilename);
        return
    end
    if isscalar(trackParams{1}) && (trackParams{1} == abs(round(trackParams{1})))
        fibMaxLen= trackParams{1};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbFACT(error): fibMaxLen has wrong properties', mfilename);
        return
    end
    if isscalar(trackParams{2}) && (trackParams{2} == abs(round(trackParams{2})))
        sampleNo= trackParams{2};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbFACT(error): repNo has wrong properties', mfilename);
        return
    end
    [mapAy, errStr, opLog, dataType]= private_ProbFACT(eigVec.dataAy, eigVal.dataAy, vox, seedVc, maskAy, ...
        sampleNo, fibMaxLen, vHd);
    if isempty(mapAy), return; end
elseif strcmp(algName, 'ProbKoch')
    if length(trackParams) < 3
        errStr= sprintf('%s::local_startRWAlg::ProbKoch(error): not enough parameters', mfilename);
        return
    end
    if isscalar(trackParams{1}) && (trackParams{1} == abs(round(trackParams{1})))
        fibMaxLen= trackParams{1};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbKoch(error): fibMaxLen has wrong properties', mfilename);
        return
    end
    if isscalar(trackParams{2}) && (trackParams{2} == abs(round(trackParams{2})))
        repNo= trackParams{2};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbKoch(error): repNo has wrong properties', mfilename);
        return
    end
    if isscalar(trackParams{3}) && (trackParams{3} == abs(trackParams{3}))
        exp_alpha= trackParams{3};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbKoch(error): exp_alpha has wrong properties', mfilename);
        return
    end
    [mapAy, errStr, opLog, dataType]= private_ProbKoch(eigVec.dataAy, eigVal.dataAy, vox, seedVc, maskAy, ...
        repNo, fibMaxLen, exp_alpha, vHd);
    if isempty(mapAy), return; end
elseif strcmp(algName, 'ProbTarget')
    if length(trackParams) < 4
        errStr= sprintf('%s::local_startRWAlg::ProbTarget(error): not enough parameters', mfilename);
        return
    end
    if isscalar(trackParams{1}) && (trackParams{1} == abs(round(trackParams{1})))
        fibMaxLen= trackParams{1};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbTarget(error): fibMaxLen has wrong properties', mfilename);
        return
    end
    if isscalar(trackParams{2}) && (trackParams{2} == abs(round(trackParams{2})))
        repNo= trackParams{2};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbTarget(error): repNo has wrong properties', mfilename);
        return
    end
    if isscalar(trackParams{3}) && (trackParams{3} == abs(trackParams{3}))
        exp_alpha= trackParams{3};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbTarget(error): exp_alpha has wrong properties', mfilename);
        return
    end
    if islogical(trackParams{4}) && isequal(sizeAy, size(trackParams{4}))
        targetAy= trackParams{4};
    else
        errStr= sprintf('%s::local_startRWAlg::ProbTarget(error): exp_alpha has wrong properties', mfilename);
        return
    end
    [mapAy, errStr, opLog, dataType]= private_ProbTarget(eigVec.dataAy, eigVal.dataAy, vox, seedVc, maskAy, ...
        targetAy, repNo, fibMaxLen, exp_alpha, vHd);
    if isempty(mapAy), return; end
else
    errStr= sprintf('%s::local_startRWAlg(error): internal error', mfilename);
    return
end

%% generiere des mapNames
mapName= sprintf('%s(%s)', algName, datestr(now));
%% erzeugen der probStruct
if strcmp(dataType, 'PMap')
    [res, errStr]= probstruct_init(mrstruct_init('volume', mapAy, dtdstruct_query(dtdStruct, 'mrStructProb')), ...
        mapName, dataType, algName, opLog);
else
    [res, errStr]= probstruct_init(mrstruct_init('series3D', mapAy, dtdstruct_query(dtdStruct, 'mrStructProb')), ...
        mapName, dataType, algName, opLog);
end
    



%%%%%%%%% B= 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%55%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function verStr= private_catchAlgVersionName(dllName)
verStr= '';

defName= strcat(dllName, '_V0.1');
try %try to catsch version name of randomwalk alg
    verStr= eval(dllNameCell{i});
catch
    verStr= [];
end

if isempty(verStr) || ~ischar(verStr)% version alg seems not available .... use default name
    verStr= defName;
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, fracAy]= private_connWeightMult(data_1, data_2, B)
res= []; errStr= ''; fracAy= [];

% find max of maps
m1= max(max(max(sum(data_1, 4))));
m2= max(max(max(sum(data_2, 4))));

% calculate connected and merged fraction
par= sum(data_1 .* data_2, 4)/(m1*m2);
apar= sum(data_1 .* data_2(:, :, :, [2 1]), 4)/(m1*m2);

% gen correspons to a simple multiplication
gen= par + apar;
idx= find(gen ~= 0);
isCon= apar(idx)./(apar(idx) + par(idx));

fracAy= zeros(size(gen));
fracAy(idx)= isCon;

% calculate the empirical probability of an underlying connection 
tran= 1 - 1./(exp((isCon - 0.5)./B) + 1);

% weight simple mult by the probability of a connectioed /merged  configuration
res= zeros(size(data_1));
res(idx)= gen(idx).*tran;                           % connection weighted mult
res(numel(gen) + idx)= gen(idx).*(1 - tran);   % merged weighted mult


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, opLog, dataType]= private_ProbRand(eigVec, eigVal, vox, seedVc, maskAy, repNo, fibMaxLen, exp_alpha, vHd)
res= []; errStr= ''; opLog= []; dataType= 'PMap';
% start alg
tic;
try
    res= random_walk_prob(eigVec, eigVal, vox, double(maskAy), seedVc, repNo, fibMaxLen, exp_alpha, zeros(size(maskAy)), vHd);        %%% detect invalid entry string
catch
    errStr= sprintf('%s::private_ProbRand(error):algorithm failed', mfilename);
    return
end

opLog.type= 'rwAlg';
opLog.date= datestr(now);
opLog.params{1}.MaxFiberLen= fibMaxLen;
opLog.params{1}.repititionNo= repNo;
opLog.params{1}.alpha= exp_alpha;
opLog.params{1}.duration= toc;
opLog.params{1}.maskAy= maskAy;
opLog.params{1}.seedVc= seedVc;
opLog.params{1}.algName= private_catchAlgVersionName('random_walk_prob');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, opLog, dataType]= private_ProbRandExt(eigVec, eigVal, vox, seedVc, maskAy, repNo, fibMaxLen, exp_alpha, allowRevisit, vHd)
res= []; errStr= ''; opLog= []; dataType= '2ModPMap';
% start alg
tic;
try
    if allowRevisit
        res= random_walk_prob(eigVec, eigVal, vox, double(maskAy), seedVc, repNo, fibMaxLen, exp_alpha, zeros([size(maskAy) 2]), vHd);
    else
        res= random_walk_prob(eigVec, eigVal, vox, double(maskAy), seedVc, repNo, fibMaxLen, exp_alpha, zeros([size(maskAy) 2]), vHd, 1);
    end
catch
    errStr= sprintf('%s::private_ProbRandExt(error):algorithm failed', mfilename);
    return
end

opLog.type= 'rwAlg';
opLog.date= datestr(now);
opLog.params{1}.MaxFiberLen= fibMaxLen;
opLog.params{1}.repititionNo= repNo;
opLog.params{1}.alpha= exp_alpha;
opLog.params{1}.duration= toc;
opLog.params{1}.maskAy= maskAy;
opLog.params{1}.seedVc= seedVc;
opLog.params{1}.algName= private_catchAlgVersionName('random_walk_prob');
opLog.params{1}.AllowRevisits= allowRevisit;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, opLog, dataType]= private_ProbRandExtHARDI(eigVec, eigVal, HARDI, vox, seedVc, maskAy, repNo, fibMaxLen, exp_alpha, allowRevisit, vHd)
res= []; errStr= ''; opLog= []; dataType= '2ModPMap';
% start alg
tic;
try
    if allowRevisit
        res= random_walk_prob_hardi(eigVec, eigVal, vox, HARDI, seedVc, repNo, fibMaxLen, exp_alpha, zeros([size(maskAy) 2]), vHd);
    else
        res= random_walk_prob_hardi(eigVec, eigVal, vox, HARDI, seedVc, repNo, fibMaxLen, exp_alpha, zeros([size(maskAy) 2]), vHd, 1);
    end
catch
    errStr= sprintf('%s::private_ProbRandExt(error):algorithm failed', mfilename);
    return
end

opLog.type= 'rwAlg';
opLog.date= datestr(now);
opLog.params{1}.MaxFiberLen= fibMaxLen;
opLog.params{1}.repititionNo= repNo;
opLog.params{1}.alpha= exp_alpha;
opLog.params{1}.duration= toc;
opLog.params{1}.maskAy= maskAy;
opLog.params{1}.seedVc= seedVc;
opLog.params{1}.algName= private_catchAlgVersionName('random_walk_prob');
opLog.params{1}.AllowRevisits= allowRevisit;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, opLog, dataType]= private_ProbRandRef(eigVec, eigVal, vox, seedVc, maskAy, repNo, fibMaxLen, exp_alpha, vHd)
res= []; errStr= ''; opLog= []; dataType= 'PMap';

errStr= sprintf('%s::private_ProbRandRef(error):algorithm is not supported yet', mfilename);
return
% start alg
tic;
try
    res= random_walk_ref(eigVec, eigVal, vox, maskAy, seedVc, repNo, fibreMaxLen, exp_alpha, zeros(size(maskAy)), vHd);        %%% detect invalid entry string
catch
    errStr= sprintf('%s::private_ProbRandRef(error):algorithm failed', mfilename);
    return
end

opLog.type= 'rwAlg';
opLog.date= datestr(now);
opLog.params{1}.MaxFiberLen= fibMaxLen;
opLog.params{1}.repititionNo= repNo;
opLog.params{1}.alpha= exp_alpha;
opLog.params{1}.duration= toc;
opLog.params{1}.maskAy= maskAy;
opLog.params{1}.seedVc= seedVc;
opLog.params{1}.algName= private_catchAlgVersionName('random_walk_ref');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, opLog, dataType]= private_ProbFACT(eigVec, eigVal, vox, seedVc, maskAy, sampleNo, fibMaxLen, vHd)
res= []; errStr= ''; opLog= []; dataType= 'PMap';

errStr= sprintf('%s::private_ProbFACT(error):algorithm is not supported yet', mfilename);
return

% start alg
tic;
try
    res= prob_fakt_main(eigVec, eigVal, vox, maskAy, seedVc, NaN, fibMaxLen, sampleNo, zeros(size(maskAy)), vHd);
catch
    errStr= sprintf('%s::private_ProbFACT(error):algorithm failed', mfilename);
    return
end

opLog.type= 'rwAlg';
opLog.date= datestr(now);
opLog.params{1}.MaxFiberLen= fibMaxLen;
opLog.params{1}.sampleNo= sampleNo;
opLog.params{1}.duration= toc;
opLog.params{1}.maskAy= maskAy;
opLog.params{1}.seedVc= seedVc;
opLog.params{1}.algName= private_catchAlgVersionName('prob_fakt_main');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, opLog, dataType]= private_ProbTarget(eigVec, eigVal, vox, seedVc, maskAy, targetAy, repNo, fibMaxLen, exp_alpha, vHd)
res= []; errStr= ''; opLog= []; dataType= 'PMap';

errStr= sprintf('%s::private_ProbTarget(error):algorithm is not supported yet', mfilename);
return

% start alg
tic;
[frak, dirVc]= local_getdir_fraction_cube(vox);
try
    res= random_walk_goal(eigVec, eigVal, vox, double(maskAy)+double(targetAy), seedVc, repNo, fibMaxLen, exp_alpha, zeros(sizeAy), vHd);        
catch
    errStr= sprintf('%s::private_ProbTarget(error):algorithm failed', mfilename);
    return
end

opLog.type= 'rwAlg';
opLog.date= datestr(now);
opLog.params{1}.MaxFiberLen= fibMaxLen;
opLog.params{1}.repititionNo= repNo;
opLog.params{1}.alpha= exp_alpha;
opLog.params{1}.duration= toc;
opLog.params{1}.maskAy= maskAy;
opLog.params{1}.targetAy= targetAy;
opLog.params{1}.seedVc= seedVc;
opLog.params{1}.algName= private_catchAlgVersionName('random_walk_goal');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, opLog, dataType]= private_ProbKoch(eigVec, eigVal, vox, seedVc, maskAy, repNo, fibMaxLen, exp_alpha, vHd)
res= []; errStr= ''; opLog= []; dataType= 'PMap';

errStr= sprintf('%s::private_ProbKoch(error):algorithm is not supported yet', mfilename);
return

% start alg
tic;
[frak, dirVc]= local_getdir_fraction_cube(vox);
try
    res= koch_alg_c(eigVec, eigVal, maskAy, seedVc, repNo, fibMaxLen, exp_alpha, dirVc, frak, zeros(size(maskAy)), vHd);
catch
    errStr= sprintf('%s::private_ProbKoch(error):algorithm failed', mfilename);
    return
end

opLog.type= 'rwAlg';
opLog.date= datestr(now);
opLog.params{1}.MaxFiberLen= fibMaxLen;
opLog.params{1}.repititionNo= repNo;
opLog.params{1}.alpha= exp_alpha;
opLog.params{1}.duration= toc;
opLog.params{1}.maskAy= maskAy;
opLog.params{1}.seedVc= seedVc;
opLog.params{1}.algName= private_catchAlgVersionName('koch_alg_c');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, dirVc]= local_getdir_fraction_cube(vox)
   
dirVc= [[1 1 1]; ...  %01
        [1 1 0]; ...  %02 
        [1 1 -1]; ... %03 
        [1 0 1]; ...  %04
        [1 0 0]; ...  %05
        [1 0 -1]; ... %06
        [1 -1 1]; ... %07
        [1 -1 0]; ... %08
        [1 -1 -1]; ...%09
        [0 1 1]; ...  %10
        [0 1 0]; ...  %11
        [0 1 -1]; ... %12
        [0 0 1]; ...  %13
        [0 0 -1]; ... %14
        [0 -1 1]; ... %15
        [0 -1 0]; ... %16
        [0 -1 -1]; ...%17
        [-1 1 1]; ... %18
        [-1 1 0]; ... %19
        [-1 1 -1]; ...%20
        [-1 0 1]; ... %21
        [-1 0 0]; ... %22
        [-1 0 -1]; ...%23
        [-1 -1 1]; ...%24
        [-1 -1 0]; ...%25
        [-1 -1 -1]];  %26
if length(vox) == 4
    vox(3)= vox(3)+vox(4);
end

len= sqrt(sum(dirVc.^2, 2));
dirVc= dirVc.*(ones(size(dirVc,1), 1)*reshape(vox(1:3), [1 3]))./(len*ones(1, 3));
len= sqrt(sum(dirVc.^2, 2));
dirVc= dirVc./(len*ones(1, 3));
res= ones(1, 26);
% ?bertragen von achtel auf ganze Kugel
% res([1 3 7 9 18 20 24 26])= .25*vox(1)*vox(2) + .25*vox(1)*vox(3) + .25*vox(2)*vox(3);
% res([10 12 15 17])= .5*vox(1)*(vox(2) + vox(3));
% res([4 6 21 23])= .5*vox(2)*(vox(1) + vox(3));
% res([2 8 19 25])= .5*vox(3)*(vox(1) + vox(2));
% res([5 22])= vox(2)*vox(3);
% res([11 16])= vox(1)*vox(3);
% res([13 14])= vox(1)*vox(2);
% res= res/sum(res);