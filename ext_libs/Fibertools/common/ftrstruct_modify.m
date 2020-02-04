function [res, errStr, res2]= ftrstruct_modify(varargin)
%
%   function [res, errStr, res2]= ftrstruct_modify(ftrStruct, command, parm1[, parm2[, ...[, parmN]]])
%
%   ftrStruct: should contain a valid ftrStruct
%   command:   spezifies the action and should contain a string like
%       {'setVox' | 'setPatient' | 'setDTD_Type' | 'insertFiber' | 'renameFiber' | 'selectFiber' | 'removeFiber'}
%   param[1..end] the parameter specified by the command string
%
% return values:
%   return: the modified ftrStruct or if an error occured the old one
%   errStr: a message string, specifying the error which occured
%   res2: additional output argument, depending on the command
% 
% Description of the different commands: 
% 'setVox'        set the voxel size
%                 param1: [1x4] array containing the voxel size [voxX voxY voxZ gapZ]
% 'setPatient'    set the patient name in the ftrStruct
%                 param1: string containging the identifier of the patient or subject
% 'insertFiber'   method to insert a fiber subset into the ftrstruct
%                 param1: fiberStruct containing the ids of the tracks, and a fibersubset name
%                 param2: string if containing 'force', subsets with the same name will be overwritten
%                 res2: string containing the identifier of the inserted subset (fiberStruct)
% 'removeFiber'   method to remove a fiber subset from the ftrstruct
%                 param1: string defining the fiber subset to delete
%                 param2: string if containing 'force', subsets with the same name will be overwritten
% 'renameFiber'   method to rename a fiber subset (fiberStruct)
%                 param1: string defining the old identifyer
%                 param2: string defining the new subset identifier
%                 res2: string containing the identifier of the inserted subset (fiberStruct)
% 'selectFiber'   method to define subset of fibers touching a mask. The new fiber subset will be inserted into the ftrStruct.
%                 Comment: this method uses the 'inArea' command of ftrstruct_query
%                 param1: mrStruct, containing the mask to select the fiber subset
%                 param2: string to identify the starting fro fiber subset
%                 param3: string , defining the identifier of the new subset
%                 param4: string if containing 'force', subsets with the same name will be overwritten
%
%  Bjoern W. Kreher
%  08/02
%
%  UNIX
%

res= [];    res2= [];   errStr= '';

maxArgNo= 10;

if length(varargin) < 2
    errStr= 'ftrstruct_modify: There have to be at least two parameters';
    return;
end

if (length(varargin) > 0) & ftrstruct_istype(varargin{1})
    ftrStruct= varargin{1};
else
    errStr= 'ftrstruct_modify: First param have to be ftrStruct';
    return;
end

if ~isstr(varargin{2})
    errStr= 'ftrstruct_modify: Command have to be a string';
    return;
else
    commandStr= varargin{2};
end

param= cell(maxArgNo, 1);
for i= 1:maxArgNo
    if length(varargin) < (i + 2)
        param{i}= [];
    else
        param{i}= varargin{i + 2};
    end
end
res= ftrStruct;
if strcmp(commandStr, 'setVox')
    [res, errStr]= local_setVox(ftrStruct, param{1});
elseif strcmp(commandStr, 'setPatient')
    [res, errStr]= local_setPatient(ftrStruct, param{1});
elseif strcmp(commandStr, 'setDTD_Type')
    [res, errStr]= local_setDTD_Type(ftrStruct, param{1});
elseif strcmp(commandStr, 'insertFiber')
    [res, errStr, res2]= local_insertFiber(ftrStruct, param{1}, param{2}, param{3});
elseif strcmp(commandStr, 'renameFiber')
    [res, errStr, res2]= local_renameFiber(ftrStruct, param{1}, param{2});
elseif strcmp(commandStr, 'selectFiber')
    [res, errStr]= local_selectFiber(ftrStruct, param{1}, param{2}, param{3});
elseif strcmp(commandStr, 'removeFiber')
    [res, errStr]= local_removeFiber(ftrStruct, param{1}, param{2});
elseif strcmp(commandStr, 'repair')
    [res, errStr]= local_repairFiber(ftrStruct);
elseif strcmp(commandStr, 'clearFiber')
    [res, errStr]= local_clearFiber(ftrStruct, param{1});
elseif strcmp(commandStr, 'prepare')
    [res, errStr]= local_prepareFiber(ftrStruct);
elseif strcmp(commandStr, 'translate')
    [res, errStr]= local_translateFiber(ftrStruct, param{1});
elseif strcmp(commandStr, 'transform')
    [res, errStr]= local_transformFiber(ftrStruct, param{1});
elseif strcmp(commandStr, 'coregister')
    [res, errStr]= local_coregister(ftrStruct, param{1});
elseif strcmp(commandStr, 'hMatrix')
    [res, errStr]= local_sethMatrix(ftrStruct, param{1});
elseif strcmp(commandStr, 'moveUp')
    [res, errStr, res2]= local_moveUp(ftrStruct, param{1});
elseif strcmp(commandStr, 'moveDown')
    [res, errStr, res2]= local_moveDown(ftrStruct, param{1});
else
    errStr= sprintf('dtd_modify: Command ''%s'' is not implemented', commandStr);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, ok]= local_moveUp(ftrStruct, id)
res= ftrStruct; errStr= ''; ok= false;

[id, errStr]= private_name2id(ftrStruct, id);

if isempty(id) | (id <= 1)
    return
end

tmp= ftrStruct.fiber{id};
ftrStruct.fiber{id}= ftrStruct.fiber{id - 1};
ftrStruct.fiber{id - 1}= tmp;

ok= true;
res= ftrStruct;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, ok]= local_moveDown(ftrStruct, id)
res= ftrStruct; errStr= ''; ok= false;

[id, errStr]= private_name2id(ftrStruct, id);

if isempty(id) | (id >= length(ftrStruct.fiber))
    return
end

tmp= ftrStruct.fiber{id};
ftrStruct.fiber{id}= ftrStruct.fiber{id + 1};
ftrStruct.fiber{id + 1}= tmp;

ok= true;
res= ftrStruct;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_sethMatrix(ftrStruct, edges)
res= ftrStruct; errStr= '';

if ~isequal(size(edges), [4 4])
    errStr= sprintf('%s::setHMatrix(error): no valid edges', mfilename);
    return
end

res.hMatrix= edges;



%%
%
%% Start:
%       [res, errStr]= local_setVox(ftrStruct, vox)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ftrStruct, errStr]= local_setVox(ftrStruct, vox)
errStr= '';
ftrStruct.vox= vox;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_setVox
%


%
%
%% Start:
%       [res, errStr]= local_setVox(ftrStruct, vox)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ftrStruct, errStr]= local_setPatient(ftrStruct, patient)
errStr= '';
ftrStruct.patient= patient;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_setVox
%


%
%
%% Start:
%       [res, errStr]= local_setDTD_Type(ftrStruct, vox)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ftrStruct, errStr]= local_setDTD_Type(ftrStruct, dtdType)
errStr= '';
ftrStruct.dtdType= dtdType;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_setDTD_Type
%




%
%
%% Start:
%       [res, errStr]= local_clearFiber(ftrStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ftrStruct, errStr]= local_clearFiber(ftrStruct, force)
errStr= '';

if strcmp(force, 'force')
    buttonName= 'Yes, sure';
else
    buttonName=questdlg('Are you sure to remove all fiber subsets?', ...
        'Removing all subsets of ftrStruct', 'Yes, sure','No','No');
end

if strcmp(buttonName, 'Yes, sure');
    ftrStruct.fiber= {};
else
    errStr='fiberstruct_modify (removeFiber): aborted by user';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_clearFiber
%


%
%
%% Start:
%       [res, errStr]= local_insertFiber(ftrStruct, id)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ftrStruct, errStr, nameStr]= local_insertFiber(ftrStruct, fiberStruct, force, fibName)
errStr= ''; nameStr= '';

[ok, errStr]= ftrstruct_istype(fiberStruct, 'fiberStruct');
if ~ok
    return
end

if isempty(fibName)
    nameStr= fiberStruct.name;
else
    nameStr= fibName;
end

nameCell= ftrstruct_query(ftrStruct, 'fiberNames');
    
if isempty(nameStr)
    [nameStr, errStr]= get_unique_str(nameCell, 'fiber', 1, 0, 1);
else
    id= private_name2id(ftrStruct, nameStr);
    if ~isempty(id)
        [ftrStruct, errStr]= local_removeFiber(ftrStruct, nameStr, force);
    end
    [nameStr, errStr]= get_unique_str(nameCell, nameStr, 0, 0, 0);
end

if isempty(nameStr) | ~isempty(errStr)
    return
end

fiberStruct.name= nameStr;
ftrStruct.fiber{end+1, 1}= fiberStruct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_insertFiber
%

%
%
%% Start:
%       [res, errStr]= local_renameFiber(ftrStruct, nameOld, nameNew)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ftrStruct, errStr, nameNew]= local_renameFiber(ftrStruct, nameOld, nameNew)
errStr= '';

[id, errStr]= private_name2id(ftrStruct, nameOld);

if isempty(id)
    return;
end

nameCell= ftrstruct_query(ftrStruct, 'fiberNames');
if isempty(nameNew)
    [nameNew, errStr]= get_unique_str(nameCell, 'fiber', 1, 0, 1);
else
    [nameNew, errStr]= get_unique_str(nameCell, newName, 0, 0, 0);
end

if ~isempty(nameNew)
    ftrStruct.fiber{id, 1}.name= nameNew;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_renameFiber
%

%
%
%% Start:
%       [res, errStr]= local_removeFiber(ftrStruct, id)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ftrStruct, errStr]= local_removeFiber(ftrStruct, id, force)
errStr= '';

id= private_name2id(ftrStruct, id);

fiberNo= prod(size((ftrStruct.fiber)));
if ~isempty(id) & (id >= 1) & (id <= fiberNo)
    if strcmp(force, 'force')
        buttonName= 'Yes, sure';
    else
        buttonName=questdlg(strcat('Are you sure to remove fiber ''', ftrStruct.fiber{id}.name, ''' ?'), ...
            'Removing subset of ftrStruct', 'Yes, sure','No','No');
    end
    if strcmp(buttonName, 'Yes, sure');
        if fiberNo == 1
            ftrStruct.fiber= {};
        else
            ftrStruct.fiber(id)= [];
        end
    else
        errStr='fiberstruct_modify (removeFiber): aborted by user';
    end
else
    errStr= 'fiberstruct_modify (removeFiber): invalid id';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_removeFiber
%

%
%
%% Start:
%       [res, errStr]= local_selectFiber(ftrStruct, area, nameRef, nameNew, silent)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ftrStruct, errStr]= local_selectFiber(ftrStruct, area, nameRef, nameNew, silent)

errStr= '';
[fiberStruct, errStr]= ftrstruct_query(ftrStruct, 'inArea', nameRef);
fiberStruct.name= nameNew;
[ftrStruct, errStr]= local_insertFiber(ftrStruct, fiberStruct, silent);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_selectFiber
%


%
%
%% Start:
%       [res, errStr]= local_translateFiber(ftrStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_translateFiber(ftrStruct, transVc)

res= ftrStruct; errStr= '';


for i= 1:length(res.curveSegCell)
    res.curveSegCell{i}(:, 1)= res.curveSegCell{i}(:, 1) + transVc(1);
    res.curveSegCell{i}(:, 2)= res.curveSegCell{i}(:, 2) + transVc(2);
    res.curveSegCell{i}(:, 3)= res.curveSegCell{i}(:, 3) + transVc(3);
end

if isfield(res.user, 'orgCurveCell')
    for i= 1:length(res.user.orgCurveCell)
        res.user.orgCurveCell{i}(:, 1)= res.user.orgCurveCell{i}(:, 1) + transVc(1);
        res.user.orgCurveCell{i}(:, 2)= res.user.orgCurveCell{i}(:, 2) + transVc(2);
        res.user.orgCurveCell{i}(:, 3)= res.user.orgCurveCell{i}(:, 3) + transVc(3);
    end
end
%
%
%% Start:
%       [res, errStr]= local_transformFiber(ftrStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [res, errStr]= local_transformFiber(ftrStruct, transMx)

res= ftrStruct; errStr= '';


for i= 1:length(res.curveSegCell)
    res.curveSegCell{i}= (transMx(1:3, 1:3)*(res.curveSegCell{i}' - 1) + transMx(1:3, 4)*ones(1, size(res.curveSegCell{i}, 1)))';
end

if isfield(res.user, 'orgCurveCell')
    for i= 1:length(res.user.orgCurveCell)
        res.user.orgCurveCell{i}= (transMx(1:3, 1:3)*(res.user.orgCurveCell{i}' - 1) + transMx(1:3, 4)*ones(1, size(res.user.orgCurveCell{i}, 1)))';
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_coregister(ftrStruct, dtdStruct)
res= ftrStruct;  errStr= '';

if ~dtdstruct_istype(dtdStruct)
    errStr= sprintf('%s::coregister(error): argument is not a dtd/mrStruct', mfilename);
    return
end

hMatrix= dtdstruct_query(dtdStruct, 'hMatrix');

% check matrices
if ~isequal(size(hMatrix), [4 4]) | ~isequal(size(ftrStruct.hMatrix), [4 4])
    errStr= sprintf('%s::coregister(error): the transformation is not well defined', mfilename);
    return;
end

% restore old curves
if isfield(res.user, 'orgCurveCell')
    ftrStruct.curveSegCell= ftrStruct.user.orgCurveCell;
    ftrStruct.user= rmfield(ftrStruct.user, 'orgCurveCell');
end
% calc transformation
crMy= eye(4); crMy(1:3, 4)= -1;
transMy= inv(hMatrix*crMy)*ftrStruct.hMatrix*crMy;
% apply transformation
for i= 1:length(ftrStruct.curveSegCell)
    % first rotation
    ftrStruct.curveSegCell{i}= ftrStruct.curveSegCell{i}*(transMy(1:3, 1:3)');
    % then translate
    for k= 1:3
        ftrStruct.curveSegCell{i}(:, k)= ftrStruct.curveSegCell{i}(:, k) + transMy(k, 4);
    end
end

% set new spatial information
ftrStruct.hMatrix= hMatrix;
ftrStruct.vox= dtdstruct_query(dtdStruct, 'vox');

res= ftrStruct;

%% Start:
%       [res, errStr]= local_repairFiber(ftrStruct);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_repairFiber(ftrStruct)

res= ftrStruct; errStr= '';

orgCrv= ftrStruct.curveSegCell;
res.curveSegCell= {};
res.connectCell= {};
res.fiber= {};
lenAy= zeros(length(orgCrv), 1);
for i= 1:length(orgCrv)
    lenAy(i)= size(orgCrv{i}, 1);        
end

[dd, idx]= sort(lenAy);
outAy= zeros(length(orgCrv), 1);
for i= length(orgCrv):-1:1
    if outAy(idx(i)) == 0
        res.curveSegCell{end + 1}= orgCrv{idx(i)};
        res.connectCell{end + 1}= length(res.curveSegCell);
        cCrv= orgCrv{idx(i)};
        aoutAy(idx(i))= 1;
        for k= (i - 1):-1:1
            if outAy(idx(k)) == 0
                cmp= zeros(size(cCrv, 1), size(orgCrv{idx(k)}, 1))';
                for ii= 1:3
                    [X, Y]= meshgrid(cCrv(:, ii), orgCrv{idx(k)}(:, ii));
                    cmp= cmp + double(X == Y);
                end
                if length(find(cmp == 3)) == size(orgCrv{idx(k)}, 1)
                    outAy(idx(k))= 1;
                elseif length(find(cmp == 3)) == 0
                else
                    disp('here was somthing komisch');
                end
                    
                
            end
        end
        disp(sprintf('Repair ftrStruct %d/%d (%d -> %d)', ...
            length(orgCrv) - i + 1, length(orgCrv), length(orgCrv) - i + 1, length(res.curveSegCell)));
    end    
end


%
%
%
%% Start:
%       [res, errStr]= local_convertFiber(ftr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_prepareFiber(ftr)
errStr= ''; res= ftr;

if ~strncmp(ftr.algoName, 'GibbsTrack', length('GibbsTrack')) ...
        & ~strncmp(ftr.algoName, 'pointTrack', length('pointTrack'))
    return
end

if isfield(ftr.user, 'orgCurveCell')
    return;
end

ftr.user.orgCurveCell= ftr.curveSegCell;
crvCell= ftr.curveSegCell;

for i= 1:length(crvCell)
    stepsAy= abs(round(crvCell{i}(2:end, :)) - round(crvCell{i}(1:(end - 1), :)));
    maxStep= max(stepsAy, [], 1);
    coef= zeros(size(crvCell{i}, 1) - 1, sum(maxStep));
    dirVc= crvCell{i}(2:end, :) - crvCell{i}(1:(end - 1), :);
    p= 0;
    for k= 1:3      % ?ber alle drei dimensionen
        
        for j= 1:maxStep(k) % so viele voxel schnittpunkte es maximal gibt
            p= p + 1;
            kV= sign(dirVc(:, k))*(j - 0.5) + round(crvCell{i}(1:(end - 1), k)) - crvCell{i}(1:(end - 1), k);
            coef(:, p)= kV./dirVc(:, k);
            
        
        end
    end
    
    coef= sort(coef, 2);
    
    % f?ge zwischen punkte ein
    crvNew= zeros(size(crvCell{i}, 1) + sum(sum(stepsAy)), 3);
    p= 0;
    for j= 1:(size(crvCell{i}, 1) - 1)
        p= p + 1;
        crvNew(p, :)= crvCell{i}(j, :);
        idx= find(abs(coef(j, :)) < 1);
        if ~isempty(idx)
            idxLen= length(idx);
            crvNew(p + (1:idxLen), :)= ones(idxLen, 1)*crvCell{i}(j, :) + (ones(idxLen, 1)*dirVc(j, :)).*(coef(j, idx)'*ones(1, 3));
            p= p + length(idx);
        end
    end
    crvNew(p+1, :)= crvCell{i}(end, :);    
    
    ftr.curveSegCell{i}= crvNew(1:(p+1), :);
end

res= ftr;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%
%
%% Start:
%       [id, errStr]= private_name2id(ftrStruct, nameStr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [id, errStr]= private_name2id(ftrStruct, nameStr)
errStr= '';

if isempty(nameStr) | isnumeric(nameStr)
    id= nameStr;
    if (id < 1) | (id > length(ftrStruct.fiber))
        id= [];
        errStr= 'ftrstruct_modify (private_name2id): invalide fiber id';
    end
    return
end

id= [];
for i= 1:length(ftrStruct.fiber)
    if strcmp(ftrStruct.fiber{i}.name, nameStr)
        id= i;
        return
    end
end

if isempty(id)
    errStr= 'ftrstruct_modify (private_name2id): invalide fiber name';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: private_name2id
%
