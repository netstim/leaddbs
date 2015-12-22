function [res, errStr, res2, res3]= maskstruct_modify(varargin)
%   function [res, errStr, res2, res3]= maskstruct_modify(maskStruct, command, parm1[, parm2[, ...[, parmN]]])
%
%   maskStruct: should contain a valid maskStruct Ver2
%   command:   spezifies the action and should contain a string like
%       {'insertMask' | 'renameFiber' | 'selectFiber' | 'removeFiber'}
%   param[1..end] the parameter specified by the command string
%
% 
% 'mrStructProb'  mrStruct (mrStruct)
% 'createMask'        maskName (string)
% 'appendMaskStruct'  maskStruct2 (maskStruct)
% 'setMask'           mask (mrStruct|dataAy), maskName (string|int)
% 'setVolumeSize'     sizeAy (int[3])
% 'moveUp'            maskName (string|int)
% 'moveDown'          maskName (string|int)
% 'removeMask'        maskName (string|int), [forceFlag]
% 'renameMask'        maskName (string|int), [newName] (string)
% 'setSection'        maskName (string|int), [2x3 double]
% 'unsetSection'        maskName (string|int), [2x3 double], flag, nameStr
% 'AND'               mask1Str (string|int), mask2Str<string|int>, [newName] (string)
% 'flipLR'            mask1Str (string|int)
% 'inv'               maskStr (string|int),  [newName](string)
% 'cut'               bBox (double[2,3])
%
% return values:
%   return: the modified maskStruct or if an error occured the old one
%   errStr: a message string, specifying the error which occured
%   res2: additional output argument, depending on the command
% 
% Description of the different commands: 
% 'setVox'        set the voxel size
%                 param1: [1x4] array containing the voxel size [voxX voxY voxZ gapZ]
%
%  Bjoern W. Kreher
%  08/02
%
%  Volkmar Glauche
%  2008/08 Fix trivial mlint syntax warnings.
%
%  UNIX
%

%#ok<*NASGU>

res= [];    res2= [];   errStr= '';
%% check input param
if length(varargin) < 2
    errStr= sprintf('%s(error): There have to be at least two parameters', mfilename);
    return;
end

if ~isempty(varargin) && maskstruct_istype(varargin{1})
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

res= maskStruct;

%% split to different commands
if strcmp(commandStr, 'mrStructProb')
    if mrstruct_istype(param{1})
        mrStruct= param{1};
    else
        errStr= sprintf('%s(setMR-Properties): param have to be a mrStruct', mfilename);
        return
    end
    [res, errStr]= local_setMRProb(maskStruct, mrStruct);
elseif strcmp(commandStr, 'createMask')
    if ischar(param{1}) || isempty(param{1})
        maskName= param{1};
    else
        errStr= sprintf('%s(createMask): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_createMask(maskStruct, maskName);
elseif strcmp(commandStr, 'appendMaskStruct')
    if maskstruct_istype(param{1})
        maskStruct2= param{1};
    else
        errStr= sprintf('%s(appendMaskStruct): invalid data', mfilename);
        return
    end    
    [res, errStr]= local_appendMaskStruct(maskStruct, maskStruct2);
elseif strcmp(commandStr, 'setMask')
    if (isempty(param{1}) || isnumeric(param{1}) || islogical(param{1}) || mrstruct_istype(param{1})) && (ischar(param{2}) || isscalar(param{2}))
        mask= param{1}; maskName= param{2};
    else
        errStr= sprintf('%s(setMask): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_setMask(maskStruct, mask, maskName);
elseif strcmp(commandStr, 'setVolumeSize')
    if isnumeric(param{1})
        sizeAy= param{1};
    else
        errStr= sprintf('%s(setVolumeSize): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_setVolSize(maskStruct, sizeAy);
elseif strcmp(commandStr, 'moveUp')
    if ischar(param{1}) || isscalar(param{1})
        maskName= param{1};
    else
        errStr= sprintf('%s(moveUp): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_moveUp(maskStruct, maskName);
elseif strcmp(commandStr, 'moveDown')
    if ischar(param{1}) || isscalar(param{1})
        maskName= param{1};
    else
        errStr= sprintf('%s(moveDown): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_moveDown(maskStruct, maskName);
elseif strcmp(commandStr, 'removeMask')
    if (ischar(param{1}) || isscalar(param{1})) && (isempty(param{2}) || ischar(param{2}))
        maskName= param{1}; forceFlag= param{2};
    else
        errStr= sprintf('%s(removeMask): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_removeMask(maskStruct, maskName, forceFlag);
elseif strcmp(commandStr, 'renameMask')
    if (ischar(param{1}) || isscalar(param{1})) && (isempty(param{2}) || ischar(param{2}))
        maskName= param{1}; newName= param{2};
    else
        errStr= sprintf('%s(renameMask): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_renameMask(maskStruct, maskName, newName);
elseif strcmp(commandStr, 'setSection')
    if (ischar(param{1}) || isscalar(param{1})) && (isnumeric(param{2}) && isequal(size(param{2}), [2 3])) ...
            && (isempty(param{3}) || ischar(param{3}))
        maskName= param{1}; secVc= param{2}; newMaskStr= param{3};
    else
        errStr= sprintf('%s(setSection): invalid data', mfilename);
        return
    end    
    [res, errStr, res2, res3]= local_setSection(maskStruct, maskName, secVc, true, newMaskStr);
elseif strcmp(commandStr, 'unsetSection')
    if (ischar(param{1}) || isscalar(param{1})) && (isnumeric(param{2}) && isequal(size(param{2}), [2 3])) ...
            && (isempty(param{3}) || ischar(param{3}))
        maskName= param{1}; secVc= param{2}; newMaskStr= param{3};
    else
        errStr= sprintf('%s(unserSection): invalid data', mfilename);
        return
    end    
    [res, errStr, res2, res3]= local_setSection(maskStruct, maskName, secVc, false, newMaskStr);
elseif strcmp(commandStr, 'AND')
    if (ischar(param{1}) || isscalar(param{1})) && (ischar(param{2}) || isscalar(param{2})) && (isempty(param{3}) || ischar(param{3}))
        mask1Str= param{1}; mask2Str= param{2}; newName= param{3};
    else
        errStr= sprintf('%s(renameMask): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_OP2(maskStruct, mask1Str, mask2Str, newName, 'AND');
elseif strcmp(commandStr, 'OR')
    if (ischar(param{1}) || isscalar(param{1})) && (ischar(param{2}) || isscalar(param{2})) && (isempty(param{3}) || ischar(param{3}))
        mask1Str= param{1}; mask2Str= param{2}; newName= param{3};
    else
        errStr= sprintf('%s(renameMask): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_OP2(maskStruct, mask1Str, mask2Str, newName, 'OR');
elseif strcmp(commandStr, 'XOR')
    if (ischar(param{1}) || isscalar(param{1})) && (ischar(param{2}) || isscalar(param{2})) && (isempty(param{3}) || ischar(param{3}))
        mask1Str= param{1}; mask2Str= param{2}; newName= param{3};
    else
        errStr= sprintf('%s(renameMask): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_OP2(maskStruct, mask1Str, mask2Str, newName, 'XOR');
elseif strcmp(commandStr, 'inv')
    if (ischar(param{1}) || isscalar(param{1})) && (isempty(param{2}) || ischar(param{2}))
        maskStr= param{1}; newName= param{2};
    else
        errStr= sprintf('%s(renameMask): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_OP1_inv(maskStruct, maskStr, newName);
elseif strcmp(commandStr, 'flipLR')
    if (ischar(param{1}) || isscalar(param{1})) && (isempty(param{2}) || ischar(param{2}))
        maskStr= param{1}; newName= param{2};
    else
        errStr= sprintf('%s(renameMask): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_OP1_flip(maskStruct, maskStr, newName);    
elseif strcmp(commandStr, 'regionGrowing')
    if (ischar(param{1}) || isscalar(param{1})) && isnumeric(param{2}) && (numel(param{2}) == 3) ...
            && (isempty(param{3}) || ischar(param{3}))
        maskStr= param{1}; posVc= reshape(param{2}, [1 3]); newName= param{3};
    else
        errStr= sprintf('%s(regionGrowing): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_regionGrow(maskStruct, maskStr, posVc, newName);
elseif strcmp(commandStr, 'growSpheric')
    if (ischar(param{1}) || isscalar(param{1})) && (isempty(param{2}) || ischar(param{2})) ...
            && (~isempty(param{3}) && isnumeric(param{3}) && (numel(param{3}) == 1) && (abs(param{3}) == param{3})) ...
            && (isempty(param{4}) || (isnumeric(param{4}) && (numel(param{4}) == 1)))
        maskStr= param{1}; newName= param{2};   rad= param{3}; vHd= param{4};
    else
        errStr= sprintf('%s(renameMask): invalid data', mfilename);
        return
    end    
    [res, errStr, res2]= local_growSpheric(maskStruct, maskStr, newName, rad, vHd);
elseif strcmp(commandStr, 'erosion')
    if (ischar(param{1}) || isscalar(param{1})) && (isempty(param{2}) || ischar(param{2})) ...
            && (isempty(param{3}) || ischar(param{3})) ...
            && (isempty(param{4}) || (isnumeric(param{4}) && (numel(param{4}) == 1)))
        maskStr= param{1};  newName= param{2};  
        planeStr= param{3}; repNo= param{4};
    else
        errStr= sprintf('%s(renameMask): invalid data', mfilename);
        return
    end        
    [res, errStr, res2]= local_morphOP_ErosDila(maskStruct, maskStr, 'erosion', planeStr, repNo, newName);
elseif strcmp(commandStr, 'dilatation')
    if (ischar(param{1}) || isscalar(param{1})) && (isempty(param{2}) || ischar(param{2})) ...
            && (isempty(param{3}) || ischar(param{3})) ...
            && (isempty(param{4}) || (isnumeric(param{4}) && (numel(param{4}) == 1)))
        maskStr= param{1};  newName= param{2};  
        planeStr= param{3}; repNo= param{4};
    else
        errStr= sprintf('%s(renameMask): invalid data', mfilename);
        return
    end        
    [res, errStr, res2]= local_morphOP_ErosDila(maskStruct, maskStr, 'dilatation', planeStr, repNo, newName);
elseif strcmp(commandStr, 'shiftUp')
    if (ischar(param{1}) || isscalar(param{1})) && (isempty(param{2}) || ischar(param{2})) ...
            && (isempty(param{3}) || (isnumeric(param{3}) && (numel(param{3}) == 1)))
        maskStr= param{1};  newName= param{2};  
        voxNo= param{3};
    else
        errStr= sprintf('%s(renameMask): invalid data', mfilename);
        return
    end        
    [res, errStr, res2]= local_OP_shift(maskStruct, maskStr, 'shiftUp', voxNo, newName);    
 elseif strcmp(commandStr, 'shiftDown')
    if (ischar(param{1}) || isscalar(param{1})) && (isempty(param{2}) || ischar(param{2})) ...
            && (isempty(param{3}) || (isnumeric(param{3}) && (numel(param{3}) == 1)))
        maskStr= param{1};  newName= param{2};  
        voxNo= param{3};
    else
        errStr= sprintf('%s(renameMask): invalid data', mfilename);
        return
    end        
    [res, errStr, res2]= local_OP_shift(maskStruct, maskStr, 'shiftDown', voxNo, newName);     
elseif strcmp(commandStr, 'shiftLeft')
    if (ischar(param{1}) || isscalar(param{1})) && (isempty(param{2}) || ischar(param{2})) ...
            && (isempty(param{3}) || (isnumeric(param{3}) && (numel(param{3}) == 1)))
        maskStr= param{1};  newName= param{2};  
        voxNo= param{3};
    else
        errStr= sprintf('%s(renameMask): invalid data', mfilename);
        return
    end        
    [res, errStr, res2]= local_OP_shift(maskStruct, maskStr, 'shiftLeft', voxNo, newName);
elseif strcmp(commandStr, 'shiftRight')
    if (ischar(param{1}) || isscalar(param{1})) && (isempty(param{2}) || ischar(param{2})) ...
            && (isempty(param{3}) || (isnumeric(param{3}) && (numel(param{3}) == 1)))
        maskStr= param{1};  newName= param{2};  
        voxNo= param{3};
    else
        errStr= sprintf('%s(renameMask): invalid data', mfilename);
        return
    end        
    [res, errStr, res2]= local_OP_shift(maskStruct, maskStr, 'shiftRight', voxNo, newName);    
elseif strcmp(commandStr, 'cut')
    if isnumeric(param{1}) && isequal(size(param{1}), [2 3]) && isequal(round(param{1}), param{1})
        bBox= param{1};  
    else
        errStr= sprintf('%s(cut): invalid data', mfilename);
        return
    end        
    [res, errStr]= local_cutMask(maskStruct, bBox);
elseif strcmp(commandStr, 'coregister')
    if dtdstruct_istype(param{1}) && (isempty(param{2}) || ishandle(param{2}))
        dtdStruct= param{1};  vHd= param{2};
    else
        errStr= sprintf('%s(cut): invalid data', mfilename);
        return
    end        
    [res, errStr, res2]= local_coregister(maskStruct, dtdStruct, vHd);
else
    errStr= sprintf('%s(error): Command ''%s'' is not implemented', mfilename, commandStr);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, ok]= local_coregister(maskStruct, dtdStruct, vHd)
res= maskStruct; errStr= ''; ok= false;

[cn, errStr]= dtdstruct_query(dtdStruct, 'cmpDimensions', ...
    mrstruct_init('volume', true(maskStruct.sizeAy), maskStruct.mrsProp));

if cn < 0 % scheint nicht m???glich
    return
elseif cn > 0   % daten passen schon aufeinander ... nix mus gemacht werden
    ok= true;
    return
end

% do coregistration
for i= 1:length(maskStruct.maskCell)
    if ~isempty(vHd)
        set(vHd, 'String', sprintf('%s: performing coregistration on ''%s''', ...
            mfilename, maskStruct.maskNamesCell{i}));
        drawnow;
    end
    mr= mrstruct_init('volume', maskStruct.maskCell{i}, maskStruct.mrsProp);
    mr_reg= register(mr, dtdStruct, 'nearest');

    if dtdstruct_query(dtdStruct, 'cmpDimensions', mr_reg) ~= 1
        errStr= sprintf('%s::coregister(error): An error occured during coregistration', mfilename);
        return;
    end
    maskStruct.maskCell{i}= mr_reg.dataAy == 1;    
end

maskStruct.mrsProp= dtdstruct_query(dtdStruct, 'mrStructProb');
maskStruct.sizeAy= dtdstruct_query(dtdStruct, 'sizeAy');

res= maskStruct;
ok= true;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, nNameStr]= local_morphOP_ErosDila(mStruct, maskStr, morphOP, morphPlane, repNo, nNameStr)
errStr= ''; res= mStruct; 

[id, errStr]= private_name2id(mStruct, maskStr);
if isempty(id);    return;    end

%%%%% kreiere neuer name und maske
if isempty(nNameStr)
    nNameStr= maskStr;       % wenn keine neue maskenname definiert ist, ???berschreibe alte
    id2= id;
else
    [res, errStr, nNameStr, id2]= local_createMask(mStruct, nNameStr);
    if ~isempty(errStr); return; end
end
% select kernel
kernelAy= zeros(3, 3, 3);
if strcmp(morphOP, 'erosion') && strcmp(morphPlane, 'XY')
    kernelAy(:, 2, 2)= 1;  kernelAy(2, :, 2)= 1;
    cmpStr= '=='; cmpVal= 5;
elseif strcmp(morphOP, 'erosion') && strcmp(morphPlane, 'XZ')
    kernelAy(:, 2, 2)= 1;  kernelAy(2, 2, :)= 1;
    cmpStr= '=='; cmpVal= 5;
elseif strcmp(morphOP, 'erosion') && strcmp(morphPlane, 'YZ')
    kernelAy(2, :, 2)= 1;  kernelAy(2, 2, :)= 1;
    cmpStr= '=='; cmpVal= 5;
elseif strcmp(morphOP, 'erosion') && (isempty(morphPlane) || strcmp(morphPlane, 'XYZ'))
    kernelAy(:, 2, 2)= 1;  kernelAy(2, :, 2)= 1;  kernelAy(2, 2, :)= 1;
    cmpStr= '=='; cmpVal= 7;
elseif strcmp(morphOP, 'dilatation') && strcmp(morphPlane, 'XY')
    kernelAy(:, 2, 2)= 1;  kernelAy(2, :, 2)= 1;
    cmpStr= '>='; cmpVal= 1;
elseif strcmp(morphOP, 'dilatation') && strcmp(morphPlane, 'XZ')
    kernelAy(:, 2, 2)= 1;  kernelAy(2, 2, :)= 1;
    cmpStr= '>='; cmpVal= 1;
elseif strcmp(morphOP, 'dilatation') && strcmp(morphPlane, 'YZ')
    kernelAy(2, :, 2)= 1;  kernelAy(2, 2, :)= 1;
    cmpStr= '>='; cmpVal= 1;
elseif strcmp(morphOP, 'dilatation') && (isempty(morphPlane) || strcmp(morphPlane, 'XYZ'))
    kernelAy(:, 2, 2)= 1;  kernelAy(2, :, 2)= 1;  kernelAy(2, 2, :)= 1;
    cmpStr= '>='; cmpVal= 1;
else
    errStr= sprintf('%s::local_morphOP_ErosDila(error): undefined morphOP (''%s'') or plane', ...
        mfilename, morphOP);
    res= mStruct;   % undo
    return
end

if isempty(repNo)
    repNo= 1;
end

%%%% EIGENTLICHE OP
tmp= mrstruct_init('volume', res.maskCell{id}, mStruct.mrsProp);
for i= 1:repNo  
    [tmp, errStr]= morph_data(tmp, kernelAy, [2 2 2], cmpStr, cmpVal, 'no');
    if ~isempty(errStr)
        res= mStruct;
        return
    end
end

if ~isempty(tmp)
    res.maskCell{id2}= tmp.dataAy ~= 0;
else
    res= mStruct;   % undo
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, nNameStr]= local_OP_shift(mStruct, maskStr, OP, voxNo, nNameStr)
errStr= ''; res= mStruct; 

[id, errStr]= private_name2id(mStruct, maskStr);
if isempty(id);    return;    end

%%%%% kreiere neuer name und maske
if isempty(nNameStr)
    nNameStr= maskStr;       % wenn keine neue maskenname definiert ist, ???berschreibe alte
    id2= id;
else
    [res, errStr, nNameStr, id2]= local_createMask(mStruct, nNameStr);
    if ~isempty(errStr); return; end
end


%%%% EIGENTLICHE OP
tmp= mrstruct_init('volume', false(size(res.maskCell{id})), mStruct.mrsProp);
sizeX = size(res.maskCell{id},1);
sizeY = size(res.maskCell{id},2);
sizeZ = size(res.maskCell{id},3);
inds = find(res.maskCell{id}==1);
inds = reshape_index(inds,[sizeX sizeY sizeZ]);
for m = 1 : size(inds,1)
    
    if strcmp(OP, 'shiftRight')
        if inds(m,2)+voxNo <= sizeY
            tmp.dataAy(inds(m,1),inds(m,2)+voxNo,inds(m,3)) = true;
        else
            tmp.dataAy(inds(m,1),inds(m,2),inds(m,3)) = true;
        end
    elseif strcmp(OP, 'shiftLeft') % shift left
        if (inds(m,2)-voxNo) > 0
            tmp.dataAy(inds(m,1),inds(m,2)-voxNo,inds(m,3)) = true;
        else
            tmp.dataAy(inds(m,1),inds(m,2),inds(m,3)) = true;
        end
    elseif strcmp(OP, 'shiftUp') % shift left
        if (inds(m,2)-voxNo) > 0
            tmp.dataAy(inds(m,1)-voxNo,inds(m,2),inds(m,3)) = true;
        else
            tmp.dataAy(inds(m,1),inds(m,2),inds(m,3)) = true;
        end   
    elseif strcmp(OP, 'shiftDown')
        if inds(m,2)+voxNo <= sizeX
            tmp.dataAy(inds(m,1)+voxNo,inds(m,2),inds(m,3)) = true;
        else
            tmp.dataAy(inds(m,1),inds(m,2),inds(m,3)) = true;
        end
    end
    
    if ~isempty(errStr)
        res= mStruct;
        return
    end
end

if ~isempty(tmp)
    res.maskCell{id2}= tmp.dataAy ~= 0;
else
    res= mStruct;   % undo
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, nNameStr]= local_regionGrow(mStruct, maskStr, pos, nNameStr)
errStr= ''; res= mStruct; 

[id, errStr]= private_name2id(mStruct, maskStr);
if isempty(id);    return;    end

%%%%% kreiere neuer name und maske
if isempty(nNameStr)
    nNameStr= maskStr;       % wenn keine neue maskenname definiert ist, ???berschreibe alte
    id2= id;
else
    [res, errStr, nNameStr, id2]= local_createMask(mStruct, nNameStr);
    if ~isempty(errStr); return; end
end
%%%% EIGENTLICHE OP
[tmp, errStr]= region_growing(mrstruct_init('volume', res.maskCell{id}, mStruct.mrsProp), '3D', pos);

if ~isempty(tmp)
    res.maskCell{id2}= tmp.dataAy ~= 0;
else
    res= mStruct;   % undo
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, nameStr]= local_growSpheric(mStruct, maskStr, nameStr, rad, vHd)
errStr= '';
res= mStruct;

[id, errStr]= private_name2id(mStruct, maskStr);
if isempty(id);    return;    end

%%%%% kreiere neuer name und maske
if isempty(nameStr)
    nameStr= maskStr;       % wenn keine neue maskenname definiert ist, ???berschreibe alte
    id2= id;
else
    [res, errStr, nameStr, id2]= local_createMask(mStruct, nameStr);
    if ~isempty(errStr); return; end
end
%%%% EIGENTLICHE OP
[tmp, errStr]= morph_data_ext(mrstruct_init('volume', res.maskCell{id}, mStruct.mrsProp), 'growByCircle', rad, vHd);
if ~isempty(tmp)
    res.maskCell{id2}= tmp.dataAy;
else
    res= mStruct;   % undo
end
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, nameStr]= local_OP1_inv(mStruct, maskStr, nameStr)
errStr= '';
res= mStruct;

[id, errStr]= private_name2id(mStruct, maskStr);
if isempty(id);    return;    end

%%%%% kreiere neuer name und maske
if isempty(nameStr)
    nameStr= maskStr;       % wenn keine neue maskenname definiert ist, ???berschreibe alte
    id2= id;
else
    [res, errStr, nameStr, id2]= local_createMask(mStruct, nameStr);
    if ~isempty(errStr); return; end
end
%%%% EIGENTLICHE OP
res.maskCell{id2}= not(mStruct.maskCell{id});

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, nameStr]= local_OP1_flip(mStruct, maskStr, nameStr)
errStr= '';
res= mStruct;

[id, errStr]= private_name2id(mStruct, maskStr);
if isempty(id);    return;    end

%%%%% kreiere neuer name und maske
if isempty(nameStr)
    nameStr= maskStr;       % wenn keine neue maskenname definiert ist, ???berschreibe alte
    id2= id;
else
    [res, errStr, nameStr, id2]= local_createMask(mStruct, nameStr);
    if ~isempty(errStr); return; end
end
%%%% EIGENTLICHE OP
for m = 1 : size(res.maskCell{id2},3)
    res.maskCell{id2}(:,:,m)= fliplr(mStruct.maskCell{id}(:,:,m));
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, nameStr]= local_OP2(mStruct, mask1Str, mask2Str, nameStr, opStr)
errStr= '';
res= mStruct;

[id1, errStr]= private_name2id(mStruct, mask1Str);
if isempty(id1);    return;    end
[id2, errStr]= private_name2id(mStruct, mask2Str);
if isempty(id1);    return;     end

%%%%% kreiere neuer name und maske
if isempty(nameStr)
    id3= id2;
else
    [res, errStr, nameStr, id3]= local_createMask(mStruct, nameStr);
    if ~isempty(errStr); return; end
end
%%%% EIGENTLICHE OP
if strcmp(opStr, 'AND')
    res.maskCell{id3}= and(mStruct.maskCell{id1}, mStruct.maskCell{id2});
elseif strcmp(opStr, 'OR')
    res.maskCell{id3}= or(mStruct.maskCell{id1}, mStruct.maskCell{id2});
elseif strcmp(opStr, 'XOR')
    res.maskCell{id3}= xor(mStruct.maskCell{id1}, mStruct.maskCell{id2});
else
    errStr= sprintf('%s::(warning)local_OP2: Operation ''%s'' is not suported yet', ...
        mfilename, opStr);
    res= mStruct;
end
    
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mStruct, errStr, nameStr, idx]= local_setSection(mStruct, maskName, secVc, val, nameStr)
errStr= ''; idx= [];

id= [];

if ~isempty(maskName)
    [id, errStr]= private_name2id(mStruct, maskName);
    if isempty(id);    return;    end
elseif isempty(nameStr)
    nameStr= 'mask';
end
    

%%%%% kreiere gegbenenfalls neuer name und maske
if isempty(nameStr)
    id2= id;
else
    [mStruct, errStr, nameStr, id2]= local_createMask(mStruct, nameStr);
    if ~isempty(errStr); return; end
    if ~isempty(id)
        mStruct.maskCell{id2}= mStruct.maskCell{id};
    end
end

%%%%% Bestimme sektion
minVc= round(min(secVc, [], 1));   maxVc= round(max(secVc, [], 1));

idx= minVc < 1; minVc(idx)= 1;
idx= maxVc > mStruct.sizeAy; maxVc(idx)= mStruct.sizeAy(idx);

[X, Y, Z]= meshgrid(minVc(1):maxVc(1), minVc(2):maxVc(2), minVc(3):maxVc(3));
idx= reshape_index_back([reshape(X, [numel(X) 1]), reshape(Y, [numel(Y) 1]), reshape(Z, [numel(Z) 1])], mStruct.sizeAy);

%%%% eigentliche OP
mStruct.maskCell{id2}(idx)= val;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_appendMaskStruct(mStruct, mStruct2, forceFlag)
errStr= '';
res= mStruct;

[maskNo, errStr]= maskstruct_query(mStruct2, 'maskNo');
if isempty(maskNo) || maskNo == 0
    return;
end

[sizeAy, errStr]= maskstruct_query(mStruct2, 'sizeAy');
if ~isequal(mStruct.sizeAy, sizeAy) && ~isempty(mStruct.sizeAy) && ~isempty(sizeAy)
    errStr= sprintf('%s::local_appendMask(warning): maskStruct is not compatible', mfilename);
    return
end

[mNames, errStr]= maskstruct_query(mStruct2, 'maskNames');
for i= 1:maskNo
    % Erzeugen des masken names
    nameStr= mNames{i};
    while any(strcmp(res.maskNamesCell, nameStr))
        nameStr= strcat('_', nameStr);
    end
    [res, errStr, id]= maskstruct_modify(res, 'createMask', nameStr);
    if ~isempty(errStr)
        res= mStruct;
        return
    end
    [dataAy, errStr]= maskstruct_query(mStruct2, 'getMask', mNames{i}, 'noExtend');
    if ~isempty(errStr)
        res= mStruct;
        return
    end
    
    [res, errStr, id]= maskstruct_modify(res, 'setMask', dataAy, nameStr);
    if ~isempty(errStr)
        res= mStruct;
        return
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, nameStr]= local_renameMask(mStruct, maskName, nameStr)
errStr= '';
res= mStruct;

[id, errStr]= private_name2id(mStruct, maskName);
if isempty(id)
    nameStr= '';
    return;
end

if isempty(nameStr)
    [nameStr, errStr]= get_unique_str(mStruct.maskNamesCell, mStruct.maskNamesCell{id}, 0, 0, 1);
else
    [nameStr, errStr]= get_unique_str(mStruct.maskNamesCell, nameStr, 0, 0, 0);
end

if isempty(nameStr)
    return
end
res.maskNamesCell{id}= nameStr;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, nextStr]= local_removeMask(mStruct, maskName, forceFlag)
res= mStruct;   errStr= ''; nextStr= '';

[id, errStr]= private_name2id(mStruct, maskName);
if isempty(id)
    return;
end

if ~isempty(forceFlag)
    buttonName=questdlg(strcat('Are you sure to remove the mask ''', mStruct.maskNamesCell{id}, ''' ?'), ...
            'Removing mask', 'Yes, sure','No','No');
else
    buttonName= 'Yes, sure';
end

if ~strcmp(buttonName, 'Yes, sure')
    errStr= sprintf('%s::local_removeMask(warning): aborted by user', mfilename);
else
    res.maskNamesCell(id)= [];
    res.maskCell(id)= [];
    if numel(res.maskNamesCell) == 0
        res.maskNamesCell= {}; 
        res.maskCell= {}; 
        return;
    elseif numel(res.maskNamesCell) >= id
        nextStr= res.maskNamesCell{id};
    else
        nextStr= res.maskNamesCell{id-1};
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, ok]= local_moveUp(mStruct, maskName)
res= mStruct; errStr= ''; ok= 0;

[id, errStr]= private_name2id(mStruct, maskName);
if isempty(id)
    return;
elseif id == 1
    errStr= sprintf('%s::local_moveUp(warning): mask is already top of the list', mfilename);
    ok= id;
elseif (length(mStruct.maskCell) < id) || (length(mStruct.maskNamesCell) < id)
    errStr= sprintf('%s::local_moveUp(warning): internal error', mfilename);
else
    res.maskNamesCell([id - 1 id])= mStruct.maskNamesCell([id id - 1]);
    res.maskCell([id - 1 id])= mStruct.maskCell([id id - 1]);
    ok= id - 1;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, ok]= local_moveDown(mStruct, maskName)
res= mStruct; errStr= ''; ok= 0;

[id, errStr]= private_name2id(mStruct, maskName);
if isempty(id)
    return;
elseif (length(mStruct.maskCell) <= id) || (length(mStruct.maskNamesCell) <= id)
    errStr= sprintf('%s::local_moveDown(warning): mask is already last element', mfilename);
    ok= id;
else
    res.maskNamesCell([id + 1 id])= mStruct.maskNamesCell([id id + 1]);
    res.maskCell([id + 1 id])= mStruct.maskCell([id id + 1]);
    ok= id + 1;
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mStruct, errStr, ok]= local_setVolSize(mStruct, sizeAy)
errStr= ''; ok= 1;

if isempty(mStruct.sizeAy)
    mStruct.sizeAy= sizeAy;
elseif ~isequal(mStruct.sizeAy, sizeAy)
    errStr= sprintf('%s::setVolSize(error): in consistent size of volume', mfilename);
    ok= 0;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mStruct, errStr]= local_setMRProb(mStruct, mrStruct)
errStr= '';

if isempty(mrStruct.vox)
    warning('maskstruct_modify:local_setMRProb','%s::local_setMRProb(warning): vox size is not defined! Set vox to [1 1 1 0]', mfilename);
    mrStruct.vox= [1 1 1 0];
end

if isempty(mrStruct.dataAy)
    mStruct.mrsProp= mrStruct;
    return
else
    [mStruct, errStr, ok]= local_setVolSize(mStruct, mrstruct_query(mrStruct, 'sizeAy'));
    if isempty(errStr)
        mrStruct.dataAy= [];
        mStruct.mrsProp= mrStruct;
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mStruct, errStr, nameStr, id]= local_createMask(mStruct, nameStr)
errStr= ''; 

if isempty(mStruct.sizeAy)
    errStr= sprintf('%s::local_createMask(error): volume size have to be set first', mfilename);
    return
end

if isempty(nameStr)
    [nameStr, errStr]= get_unique_str(mStruct.maskNamesCell, 'mask', 1, 0, 1);
else
    [nameStr, errStr]= get_unique_str(mStruct.maskNamesCell, nameStr, 0, 0, 0);
end

if isempty(nameStr)
    return
end

id= numel(mStruct.maskNamesCell) + 1;
mStruct.maskNamesCell{id, 1}= nameStr;
mStruct.maskCell{id, 1}= false(mStruct.sizeAy);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mStruct, errStr, id]= local_setMask(mStruct, mask, maskName)
errStr= '';

[id, errStr]= private_name2id(mStruct, maskName);
if isempty(id)
    return;
end

if mrstruct_istype(mask)
    mask= mask.dataAy;
end

if isempty(mask)
    mStruct.maskCell{id}= [];
    return;
end

sizeAy= size(mask);
if isempty(mStruct.sizeAy)
    mStruct.sizeAy= sizeAy;
elseif ~isequal(mStruct.sizeAy, sizeAy)
    errStr= sprintf('%s::setMask(error): in consistent size of volume', mfilename);
    id= [];
    return
end

if islogical(mask)
    mStruct.maskCell{id}= mask;
else
    mStruct.maskCell{id}= logical(mask);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mStruct, errStr]= local_cutMask(mStruct, bBoxVc)
errStr= ''; 

bBoxVc= [min(bBoxVc, [], 1); max(bBoxVc, [], 1)];
if any(bBoxVc(1, :) < 1) || any(bBoxVc(2, :) > mStruct.sizeAy)
    errStr= sprintf('%s::local_cutMask(error): bBox exceeds the volume', mfilename);
    return
end

for i= 1:length(mStruct.maskCell)
    mStruct.maskCell{i}= mStruct.maskCell{i}(bBoxVc(1, 1):bBoxVc(2, 1), bBoxVc(1, 2):bBoxVc(2, 2), bBoxVc(1, 3):bBoxVc(2, 3));
end
mStruct.sizeAy= bBoxVc(2, :) - bBoxVc(1, :) + 1;

mStruct.user.cutingBBox= bBoxVc;
mStruct.mrsProp.user.cutingBBox= bBoxVc;
if isequal(size(mStruct.mrsProp.edges), [4, 4])
    trAy= eye(4); trAy(1:3, 4)= bBoxVc(1, :) - 1;
    mStruct.mrsProp.edges= mStruct.mrsProp.edges*trAy;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  START:
%       [id, errStr]= private_name2id(mStruct, nameStr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [id, errStr]= private_name2id(mStruct, nameStr)
id = []; 
errStr= '';

if isempty(nameStr) || isnumeric(nameStr)
    id= nameStr;
    if (id < 1) || (id > length(mStruct.maskNamesCell))
        id= [];
        errStr= 'maskstruct_modify(private_name2id): invalide mask id';
    end
    return
end

id= find(strcmp(mStruct.maskNamesCell, nameStr));

if isempty(id)
    errStr= 'maskstruct_modify(private_name2id): invalide mask name';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  END: private_name2id
