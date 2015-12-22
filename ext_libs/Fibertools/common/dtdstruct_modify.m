function [res, errStr, oArg1]= dtdstruct_modify(varargin)
%
%   function [res, errStr, oArg1]= dtdstruct_modify(dtdStruct, command, parm1[, parm2[, ...[, parmN]]])
%
%   General method to manipulate a dtdStruct
%
%   dtdStruct: should contain a valid dtdStruct
%   command:   spezifies the action and should contain a string like
%       ['setOM' | 'setVox' | 'mask' | 'user' | 'addStruct' | 'norm' | 'scale' | 'sortEigvec]
%   param[1..end] the parameter specified by the command string
%
% return values:
%   return: the modified dtdStruct or if an error occured the old one
%   errStr: a message string, specifying the error which occured
%   oArg1: additional output argument, depending on the command
% 
% Description of the different commands: 
%  'setOM' allows to set the 4x4 transformation from voxel index to world system
%         param1: homogen 4x4 matrix
%  'setVox' allows to change the vox size of the dtdStruct.
%         param1: voxel size 1x4 vector [voxX voxY voxZ gapZ]
%  'mask' multiplies the given mask with the the dtdStruct.
%         param1: a mrStruct with the same size as the dtdStruct containing the mask to apply
%  'user' allows to add a field to the user entry
%         param1: Identification string for the user entry
%         param2: data of arbitary type which will be saved under the given identification string
%  'addStruct' allows to add an additional map with of the same size to the dtdStruct
%         param1: a mrStruct of the same size as the dtdStruct, containing the map
%         param2: Identification string to identify the map
%         oArg1: the commandStr name identifying the map inside the dtdStruct
%  'norm'  Normalizes the eigenvectors to the length one
%  'scale' Scales the eigenvectors
%         param1: scaling vector
%  'sortEigvec' sorts the eigenvectors depending to the eigenvalues. 
%         By doing this, the performace of dtdstruct_query(dtd, 'getDiffDir') will be strongly increased
%
%
%  Bjoern W. Kreher
%  08/02
%
%  UNIX
%


res= varargin{1};    errStr= '';    oArg1= [];
maxArg= 10;

[isT, dtdType, errStr]= dtdstruct_istype(varargin{1});

if isT && strcmp(dtdType, 'mrStruct')
    varargin{1}= dtdstruct_init('MR', varargin{1});
end


if ischar(varargin{1}) && exist(varargin{1}, 2)
    dtdStruct= dtdstruct_read(varargin{1});
elseif isT
    dtdStruct= varargin{1};
else
    errStr= 'Error in function [res, errStr]= dtd_modify(varagrin): first param have to be dtdStruct';
    return;
end

if ~ischar(varargin{2})
    errStr= 'Error in function [res, errStr]= dtd_modify(varagrin): command have to be a string';
    return;
end
commandStr= varargin{2};

param= cell(maxArg, 1);
for i= 1:maxArg
    if length(varargin) < (i + 2)
        param{i}= [];
    else
        param{i}= varargin{i + 2};
    end
end

dtdType= dtdstruct_query(dtdStruct, 'dtdType');

if strcmp(commandStr, 'setOM')
    orM= param{1};
    [res, errStr]= local_setOM(dtdStruct, orM);
elseif strcmp(commandStr, 'setVox')
    vox= param{1};
    [res, errStr]= local_setVox(dtdStruct, vox);
elseif strcmp(commandStr, 'mask')
    mask= param{1};
    [res, errStr]= local_maskData(dtdStruct, mask);
elseif strcmp(commandStr, 'cut')
    bBox= param{1};
    [res, errStr]= local_cutData(dtdStruct, bBox);
elseif strcmp(commandStr, 'user')
    if ischar(param{1})
        userStr= param{1};
    else
        errStr= strcat(mfilename, '(error): userStr have to be of type string');
        return
    end
    dataCt= param{2};
    [res, errStr]= local_setUser(dtdStruct, userStr, dataCt);
elseif strcmp(commandStr, 'addStruct')
    mrStruct= param{1};
    structName= param{2};
    guiFlag= param{3};
    [res, errStr, oArg1]= local_addStruct(dtdStruct, structName, mrStruct, guiFlag);
elseif strcmp(commandStr, 'removeStruct')
    mapName= param{1};
    guiFlag= param{3};
    [res, errStr]= local_removeStruct(dtdStruct, mapName, guiFlag);
elseif strcmp(commandStr, 'renameStruct')
    oldName= param{1};
    newName= param{2};
    guiFlag= param{3};
    [res, errStr, oArg1]= local_renameStruct(dtdStruct, oldName, newName, guiFlag);
elseif strcmp(commandStr, 'norm')  
    [res, errStr]= local_normEigVec(dtdStruct);
elseif strcmp(commandStr, 'scale')  
    scale= param{1};
    [res, errStr]= local_scaleEigVec(dtdStruct, scale);
elseif strcmp(commandStr, 'sortEigvec')  
    if strcmp(dtdType, 'DTD')
        [dtdStruct.eigenVec_struc, dtdStruct.eigenVal_struc]= local_sortEigVec(dtdStruct.eigenVec_struc, dtdStruct.eigenVal_struc);
        res= dtdStruct;
    elseif strcmp(dtdType, 'monoMDT')
        [dtdStruct.eigenVec_struc, dtdStruct.eigenVal_struc]= local_sortEigVec(dtdStruct.eigenVec_struc, dtdStruct.eigenVal_struc);
        [dtdStruct.eigVec_1, dtdStruct.eigVal_1]= local_sortEigVec(dtdStruct.eigVec_1, dtdStruct.eigVal_1);
        res= dtdStruct;
    elseif strcmp(dtdType, 'multiDTD')
        [dtdStruct.eigenVec_struc, dtdStruct.eigenVal_struc]= local_sortEigVec(dtdStruct.eigenVec_struc, dtdStruct.eigenVal_struc);
        [dtdStruct.eigVec_1, dtdStruct.eigVal_1]= local_sortEigVec(dtdStruct.eigVec_1, dtdStruct.eigVal_1);
        [dtdStruct.eigVec_2, dtdStruct.eigVal_2]= local_sortEigVec(dtdStruct.eigVec_2, dtdStruct.eigVal_2);
        res= dtdStruct;
    else
        errStr= ['sort is not supported for dtdStruct type ''' dtdType ''''];
    end
elseif strcmp(commandStr, 'sortDT_dir') && strcmp(dtdType, 'multiDTD')
    [res, errStr]= local_sortDT_dir(dtdStruct);
elseif strcmp(commandStr, 'sortDT_x-axis') && strcmp(dtdType, 'multiDTD')
    [res, errStr]= local_sortDT_axis(dtdStruct);
elseif strcmp(commandStr, 'est_frac_1') && strcmp(dtdType, 'multiDTD')
    bValue= param{1};
    [res, errStr]= local_est_frac_1(dtdStruct, bValue);
else
    errStr= sprintf('Error in function [res, errStr]= dtd_modify(varagrin): command ''%s'' is not implemented', commandStr);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
%
%  START:
%       [res, errStr]= local_cutData(dtdStruct, bBox)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_cutData(dtdStruct, bBox)

res= dtdStruct; errStr= '';

stNames= fieldnames(dtdStruct);
dtdCell= struct2cell(dtdStruct);

sizeDTDAy= dtdstruct_query(dtdStruct, 'sizeAy');

bBox= sort(bBox);
if ~isequal(size(bBox), [2 3])
    errStr= strcat(mfilename, '::local_cutData (error): bBox have to be the size [2 3]');
    return
end
if (bBox(1, :) < 1) || (bBox(2, :) > sizeDTDAy)
    errStr= strcat(mfilename, '::local_maskData (error): bBox doesn''t fit in in dataSet');
    return;
end

sizeNewAy= 1 + bBox(2, :) - bBox(1, :);
for i= 1:length(dtdCell)
    if mrstruct_istype(dtdCell{i})
        sizCur= size(dtdCell{i}.dataAy);
        dtdCell{i}.dataAy= dtdCell{i}.dataAy(bBox(1, 1):bBox(2, 1), bBox(1, 2):bBox(2, 2), bBox(1, 3):bBox(2, 3), :);
        sizCur(1:3)= sizeNewAy;
        dtdCell{i}.dataAy= reshape(dtdCell{i}.dataAy, sizCur);
        dtdCell{i}.user.cutingBBox= bBox;
        if isequal(size(dtdCell{i}.edges), [4, 4])
            trAy= eye(4); trAy(1:3, 4)= bBoxVc(1, :) - 1;
            dtdCell{i}.edges= dtdCell{i}.edges*trAy;
        end
    end
end
res= cell2struct(dtdCell, stNames, 1);


%
%
%  START:
%       [res, errStr]= local_maskData(dtdStruct, mask)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_maskData(dtdStruct, mask)

res= dtdStruct; errStr= '';

stNames= fieldnames(dtdStruct);
dtdCell= struct2cell(dtdStruct);

sizeMaskAy= mrstruct_query(mask, 'sizeAy');
sizeDTDAy= dtdstruct_query(dtdStruct, 'sizeAy');

if ~isequal(sizeMaskAy, sizeDTDAy)
    errStr= strcat(mfilename, '::local_maskData (error): dtdStruct and mask have to be of the same size');
    return;
end

idx= find(mask.dataAy == 0);
offset= prod(sizeMaskAy);

for i= 1:length(dtdCell)
    if mrstruct_istype(dtdCell{i})
        [a, b, c, dimNo]= size(dtdCell{i}.dataAy);
    
        for j= 1:dimNo
            dtdCell{i}.dataAy(idx + (j - 1)*offset)= 0;
        end
    end
end
res= cell2struct(dtdCell, stNames, 1);
    
%
%
%  START:
%       [res, errStr]= local_setOM(dtdStruct, orM)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_setOM(dtdStruct, orM)

res= []; errStr= '';
stNames= fieldnames(dtdStruct);
dtdCell= struct2cell(dtdStruct);

for i= 1:length(dtdCell)
    dtdCell{i}.user.orientM= orM;
end

res= cell2struct(dtdCell, stNames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_setOM
%



%
%
%  START:
%       [res, errStr]= local_setUser(dtdStruct, userStr, dataCt)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_setUser(dtdStruct, userStr, dataCt)
res= dtdStruct; errStr= '';

dtdCell= struct2cell(dtdStruct);
dtdNames= fieldnames(dtdStruct);
for i= 1:length(dtdCell)
    if mrstruct_istype(dtdCell{i})
        try
            dtdCell{i}.user.(userStr) = dataCt;
        catch
            errStr= strcat(mfilename, '::local_setUser(error): An error occourred while isering dataCt in user struct');
            return
        end
    end
end

res= cell2struct(dtdCell, fieldnames(dtdStruct), 1);

%
%
%  START:
%       [res, errStr]= local_setVox(dtdStruct, vox)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_setVox(dtdStruct, vox)

res= []; errStr= '';
stNames= fieldnames(dtdStruct);
dtdCell= struct2cell(dtdStruct);

for i= 1:length(dtdCell)
    if mrstruct_istype(dtdCell{i})
        dtdCell{i}.vox= vox;
    end
end

res= cell2struct(dtdCell, stNames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_setVox
%



%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, sName]= local_renameStruct(dtdStruct, oldName, newName, guiFlag)
res= dtdStruct; errStr= ''; sName= oldName;

% remove _
if ~isempty(newName) && (newName(1) == '_')
    newName= newName(2:end);
end
if ~isempty(oldName) 
    oldName= oldName(2:end);
end

stNames= fieldnames(dtdStruct);
idx= find(strcmp(stNames, oldName));
% test if entry exist
if isempty(idx)
    errStr= sprintf('%s::renameStruct(error): the struct does not exist in the dtdStruct', mfilename);
    return;
end

% take test and delete ald structure
mrStruct= dtdStruct.(oldName);
if ~mrstruct_istype(mrStruct, 'volume')
    errStr= sprintf('%s::renameStruct(error): the struct can not be renamed', mfilename);
    return;
end
dtdStruct= rmfield(dtdStruct, oldName);
stNames{idx}= 'version'; % version is resereved

% generiere bzw zesze neuer name
if isempty(guiFlag)
    if isempty(newName)
        [newName, errStr]= get_unique_str(fieldnames(dtdStruct), oldName, 0, 1, 1);
    else
        [newName, errStr]= get_unique_str(fieldnames(dtdStruct), newName, 0, 1, 0);    
    end
    if isempty(newName)
        return;
    end
end

% Test once again of unique
if any(strcmp(stNames, newName))
    errStr= sprintf('%s::renameStruct(error): The new name does exist already', mfilename);
    return;
end
    
% weise neuer name zu
try
    dtdStruct.(newName) = mrStruct;
catch
    errStr= strcat(mfilename, '::renameStruct(error): a non valid name was specified');
	return
end
sName= strcat('_', newName);
res= dtdStruct;



%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, comStr]= local_removeStruct(dtdStruct, oldName, guiFlag)

res= dtdStruct; errStr= '';

if ~isempty(oldName) % removing the precent '_'
    oldName= oldName(2:end);
end

stNames= fieldnames(dtdStruct);
idx= find(strcmp(stNames, oldName), 1);
% test if entry exist
if isempty(idx)
    errStr= sprintf('%s::removeStruct(error): the struct does not exist in the dtdStruct', mfilename);
    return;
end

% take test and delete ald structure
mrStruct= dtdStruct.(oldName);
if ~mrstruct_istype(mrStruct, 'volume')
    errStr= sprintf('%s::removeStruct(error): the struct can not be renamed', mfilename);
    return;
end
res= rmfield(dtdStruct, oldName);



%
%
%  START:
%       [res, errStr, comStr]= local_addStruct(dtdStruct, stName, mrStruct, guiFlag)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, comStr]= local_addStruct(dtdStruct, stName, mrStruct, guiFlag)

res= dtdStruct; errStr= ''; comStr= '';
if ~mrstruct_istype(mrStruct) 
    errStr= strcat(mfilename, '::local_addStruct (error): argument was not a valid mrStruct');
    return
end

dtdSize= dtdstruct_query(dtdStruct, 'sizeAy');
mrSize= mrstruct_query(mrStruct, 'sizeAy');
if  ~isequal(dtdSize, mrSize(1:3))
    errStr= strcat(mfilename, '::local_addStruct (error): dtdStruct and mrStruct have to be of the same size');
    return
end

if isempty(guiFlag)
    if isempty(stName)
        [stName, errStr]= get_unique_str(fieldnames(dtdStruct), 'mrStruct', 1, 1, 1);
    else
        [stName, errStr]= get_unique_str(fieldnames(dtdStruct), stName, 0, 1, 0);    
    end
    if isempty(stName)
        return;
    end
    res= dtdStruct;
    res.(stName) = mrStruct;
    comStr= strcat('_', stName);
else
    nameCell= fieldnames(dtdStruct);
    try 
        res= dtdStruct;
        res.(stName) = mrStruct;
        comStr= strcat('_', stName);
    catch
        errStr= strcat(mfilename, '::local_addStruct (error): a non valid name was specified');
        res= dtdStruct;
    end
    if any(strcmp(nameCell, stName))
        res= dtdStruct;
        errStr= strcat(mfilename, '::local_addStruct (error): a non unique name was specified');
        comStr= '';
        return
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_addStruct
%

%
%
%  START:
%       function [res, errStr]= local_normEigVec(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_normEigVec(dtdStruct)

res= []; errStr= '';

stNames= fieldnames(dtdStruct);
dtdCell= struct2cell(dtdStruct);

for i= 1:length(dtdCell)
    sizeAy= size(dtdCell{i}.dataAy);
    if length(sizeAy) == 5
        vects= reshape(permute(dtdCell{i}.dataAy, [1 2 3 5 4]), [prod(sizeAy)/3 3]);
        normVals= sum(vects.^2, 2);
        idx= find(normVals ~= 0);
        normVals= 1./sqrt(normVals(idx));
        vects(idx, :)= vects(idx, :).*[normVals normVals normVals];
        dtdCell{i}.dataAy= permute(reshape(vects, sizeAy), [1 2 3 5 4]);
    end
end
res= cell2struct(dtdCell, stNames, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_normEigVec
%

%
%
%  START:
%       function [res, errStr]= local_scaleEigVec(dtdStruct, scale)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_scaleEigVec(dtdStruct, scale)

res= []; errStr= '';

stNames= fieldnames(dtdStruct);
dtdCell= struct2cell(dtdStruct);

for i= 1:length(dtdCell)
    if mrstruct_istype(dtdCell{i})
        sizeAy= size(dtdCell{i}.dataAy);
        if (length(sizeAy) == 5) && isequal(sizeAy(4:5), [3 3])
            vects= reshape(permute(dtdCell{i}.dataAy, [1 2 3 5 4]), [prod(sizeAy)/3 3]);
            vects= vects.*(ones(size(vects, 1), 1)*scale);
            normVals= sum(vects.^2, 2);
            idx= find(normVals ~= 0);
            normVals= 1./sqrt(normVals(idx));
            vects(idx, :)= vects(idx, :).*[normVals normVals normVals];
            dtdCell{i}.dataAy= permute(reshape(vects, sizeAy), [1 2 3 5 4]);
        end
    end
end
res= cell2struct(dtdCell, stNames, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_scaleEigVec
%


%
%
%  START:
%       function [res, errStr]= local_sortEigVec(eigVec, eigVal)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eigVec, eigVal]= local_sortEigVec(eigVec, eigVal)

res= []; errStr= '';

if isstruct(eigVec.user) && (sum(strcmp(fieldnames(eigVec.user), 'sort')) > 0)
    return
end

sizeAy= size(eigVec.dataAy);
voxNo= prod(sizeAy(1:3));

vects= reshape(permute(eigVec.dataAy, [4 5 1 2 3]), [3 voxNo*3]);
%vals= reshape(eigVal.dataAy, [voxNo 3]);
vals= reshape(permute(squeeze(eigVal.dataAy), [4 1 2 3]), [3 voxNo]);

[dummy, offset]= sort(vals, 1);

base= zeros(voxNo*3, 1);

base(1:3:(voxNo*3))= (0:(voxNo - 1))*3;
base(2:3:(voxNo*3))= (0:(voxNo - 1))*3;
base(3:3:(voxNo*3))= (0:(voxNo - 1))*3;

adr= base + reshape(offset, [voxNo*3 1]);

eigVec.dataAy= permute(reshape(vects(:, adr), sizeAy([4 5 1 2 3])), [3 4 5 1 2]);
eigVal.dataAy= permute(reshape(vals(adr), sizeAy([4 1 2 3])), [2 3 4 1]);

eigVec.user.sort= [3 2 1];
eigVal.user.sort= [3 2 1];


        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_sortEigVec
%



%
%
%  START:
%       function [dtdStruct, errStr]= local_sortDT_dir(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dtdStruct, errStr]= local_sortDT_dir(dtdStruct)

sizeAy= dtdstruct_query(dtdStruct, 'sizeAy');

errStr= '';
%dtdStruct.switcher = dtdStruct.va_2;
%dtdStruct.switcher.dataAy= zeros(sizeAy);

for z= 1:sizeAy(3)
    for y= 1:sizeAy(2)
        for x= 1:sizeAy(1)
            if sum(abs(dtdStruct.eigVal_1.dataAy(x, y, z, :))) ~= 0
                [dummy, maxDT_Idx]= max(squeeze(dtdStruct.eigenVal_struc.dataAy(x, y, z, :)));
                [dummy, maxDT_A_Idx]= max(squeeze(dtdStruct.eigVal_1.dataAy(x, y, z, :)));
                [dummy, maxDT_B_Idx]= max(squeeze(dtdStruct.eigVal_2.dataAy(x, y, z, :)));
                
                dirDT_Vc= squeeze(dtdStruct.eigenVec_struc.dataAy(x, y, z, :, maxDT_Idx));
                dirDT_A_Vc= squeeze(dtdStruct.eigVec_1.dataAy(x, y, z, :, maxDT_A_Idx));
                dirDT_B_Vc= squeeze(dtdStruct.eigVec_2.dataAy(x, y, z, :, maxDT_B_Idx));
                
                if abs(dirDT_Vc'*dirDT_A_Vc) < abs(dirDT_Vc'*dirDT_B_Vc)
                    tmp= squeeze(dtdStruct.eigVal_2.dataAy(x, y, z, :));
                    dtdStruct.eigVal_2.dataAy(x, y, z, :)= dtdStruct.eigVal_1.dataAy(x, y, z, :);
                    dtdStruct.eigVal_1.dataAy(x, y, z, :)= tmp;
                    
                    tmp= squeeze(dtdStruct.eigVec_2.dataAy(x, y, z, :, :));
                    dtdStruct.eigVec_2.dataAy(x, y, z, :, :)= dtdStruct.eigVec_1.dataAy(x, y, z, :, :);
                    dtdStruct.eigVec_1.dataAy(x, y, z, :, :)= tmp;
                    
%                    tmp= squeeze(dtdStruct.va_2.dataAy(x, y, z));
%                    dtdStruct.va_2.dataAy(x, y, z)= dtdStruct.va_1.dataAy(x, y, z);
%                    dtdStruct.va_1.dataAy(x, y, z)= tmp;
%                    
%                    tmp= squeeze(dtdStruct.vb_2.dataAy(x, y, z));
%                    dtdStruct.vb_2.dataAy(x, y, z)= dtdStruct.vb_1.dataAy(x, y, z);
%                    dtdStruct.vb_1.dataAy(x, y, z)= tmp;
                    
                    dtdStruct.fraction.dataAy(x, y, z, :)= dtdStruct.fraction.dataAy(x, y, z, [2 1 3]);
                    
%                    dtdStruct.switcher.dataAy(x, y, z)= rand(1, 1) + 10;
                else
%                    dtdStruct.switcher.dataAy(x, y, z)= rand(1, 1) + 1;
                end
            end
        end
    end
    z
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_sortDT_dir
%


%
%
%  START:
%       function [dtdStruct, errStr]= local_sortDT_axis(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dtdStruct, errStr]= local_sortDT_axis(dtdStruct)

sizeAy= dtdstruct_query(dtdStruct, 'sizeAy');

errStr= '';

%dtdStruct.switcher = dtdStruct.va_2;
%dtdStruct.switcher.dataAy= zeros(sizeAy);

for z= 1:sizeAy(3)
    for y= 1:sizeAy(2)
        for x= 1:sizeAy(1)
            if sum(dtdStruct.fraction.dataAy(x, y, z, :)) ~= 0
                [dummy, maxDT_Idx]= max(squeeze(dtdStruct.eigenVal_struc.dataAy(x, y, z, :)));
                [dummy, maxDT_A_Idx]= max(squeeze(dtdStruct.eigVal_1.dataAy(x, y, z, :)));
                [dummy, maxDT_B_Idx]= max(squeeze(dtdStruct.eigVal_2.dataAy(x, y, z, :)));
                
                dirDT_Vc= squeeze(dtdStruct.eigenVec_struc.dataAy(x, y, z, :, maxDT_Idx));
                dirDT_A_Vc= squeeze(dtdStruct.eigVec_1.dataAy(x, y, z, :, maxDT_A_Idx));
                dirDT_B_Vc= squeeze(dtdStruct.eigVec_2.dataAy(x, y, z, :, maxDT_B_Idx));
                
                if abs(dirDT_A_Vc(1)) > abs(dirDT_B_Vc(1))
                    tmp= squeeze(dtdStruct.eigVal_2.dataAy(x, y, z, :));
                    dtdStruct.eigVal_2.dataAy(x, y, z, :)= dtdStruct.eigVal_1.dataAy(x, y, z, :);
                    dtdStruct.eigVal_1.dataAy(x, y, z, :)= tmp;
                    
                    tmp= squeeze(dtdStruct.eigVec_2.dataAy(x, y, z, :, :));
                    dtdStruct.eigVec_2.dataAy(x, y, z, :, :)= dtdStruct.eigVec_1.dataAy(x, y, z, :, :);
                    dtdStruct.eigVec_1.dataAy(x, y, z, :, :)= tmp;
                    
                    tmp= squeeze(dtdStruct.va_2.dataAy(x, y, z));
                    dtdStruct.va_2.dataAy(x, y, z)= dtdStruct.va_1.dataAy(x, y, z);
                    dtdStruct.va_1.dataAy(x, y, z)= tmp;
                    
                    tmp= squeeze(dtdStruct.vb_2.dataAy(x, y, z));
                    dtdStruct.vb_2.dataAy(x, y, z)= dtdStruct.vb_1.dataAy(x, y, z);
                    dtdStruct.vb_1.dataAy(x, y, z)= tmp;
                    
                    dtdStruct.fraction.dataAy(x, y, z, :)= dtdStruct.fraction.dataAy(x, y, z, [2 1 3]);
                    
 %                   dtdStruct.switcher.dataAy(x, y, z)= rand(1, 1) + 10;
                else
 %                   dtdStruct.switcher.dataAy(x, y, z)= rand(1, 1) + 1;
                end
            end
        end
    end
    z
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_sortDT_axis
%

%
%
%  START:
%       function [res, errStr]= local_est_frac_1(dtdStruct);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dtdStruct, errStr]= local_est_frac_1(dtdStruct, bValue)

errStr= '';
tmp= dtdStruct.fraction.dataAy;
idx= find(tmp == 0);
tmp(idx)= 1;
dtdStruct.eigVal_1.dataAy= dtdStruct.eigVal_1.dataAy - log(cat(4, tmp(:, :, :, 1), tmp(:, :, :, 1), tmp(:, :, :, 1)))/bValue;
dtdStruct.eigVal_2.dataAy= dtdStruct.eigVal_2.dataAy - log(cat(4, tmp(:, :, :, 2), tmp(:, :, :, 2), tmp(:, :, :, 2)))/bValue;
dtdStruct.eigVal_3.dataAy= dtdStruct.eigVal_3.dataAy...
    - log(cat(4, tmp(:, :, :, 3), tmp(:, :, :, 3), tmp(:, :, :, 3)))/bValue;
dtdStruct.fraction.dataAy= zeros(size(dtdStruct.fraction.dataAy));
dtdStruct.fraction.dataAy(idx)= 1;
%dtdStruct.fraction= [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_sortDT_axis
%
