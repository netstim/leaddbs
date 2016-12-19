function [res, errStr]= dtdstruct_init(varargin)
%   function [res, errStr]= dtdstruct_init(typeStr, mrStruct, [mrData1, nameStr1, [mrData2, nameStr2, ...]])
%
%  Creates a valid dtdStruct of the type DTD
%
%typeStr 'DTD':
% [res, errStr]= dtdstruct_init(typeStr, eigenVec, eigenVal, [mrData1, nameStr1, [mrData2, nameStr2, ...]])
%   eigenVec: should contain a mrStruct of the size [sizeX sizeY sizeZ 3 3]
%      containing the eigenvectors of the diffusion tensor
%   eigenVal: should contain a mrStruct of the size [sizeX sizeY sizeZ 3]
%      containing the eigenvalues of the diffusion tensor
%   mrData[1..end]: should contain a mrStruct of the size [sizeX sizeY sizeZ]
%      it can be any arbitary map which will added to the dtdStruct
%   nameStr[1..end]: should contain a string identifiying the previous
%      mrStruct. It must be unique inside the dtdStruct and should not
%      contain special letters
%
%typeStr 'VelVd':
% [res, errStr]= dtdstruct_init(typeStr, velocityVec, [mrData1, nameStr1, [mrData2, nameStr2, ...]])
%   velocityVec: should contain a mrStruct of the size [sizeX sizeY sizeZ 3]
%      containing the eigenvalues of the diffusion tensor
%   mrData[1..end]: should contain a mrStruct of the size [sizeX sizeY sizeZ]
%      it can be any arbitary map which will added to the dtdStruct
%   nameStr[1..end]: should contain a string identifiying the previous
%      mrStruct. It must be unique inside the dtdStruct and should not
%      contain special letters
%
%
%typeStr 'MR':
% [res, errStr]= dtdstruct_init(typeStr, mrData[, nameStr], [mrData1, nameStr1, [mrData2, nameStr2, ...]])
%   mrData: should contain a mrStruct of the size [sizeX sizeY sizeZ]
%      containing the eigenvalues of the diffusion tensor
%   nameStr: should contain a string identifiying the previous
%      mrStruct. It must be unique inside the dtdStruct and should not
%      contain special letters
%   mrData[1..end]: should contain a mrStruct of the size [sizeX sizeY sizeZ]
%      it can be any arbitary map which will added to the dtdStruct
%   nameStr[1..end]: should contain a string identifiying the previous
%      mrStruct. It must be unique inside the dtdStruct and should not
%      contain special letters
%
%   return values
%     res:  contains the created dtdStruct. If an error occurs res is empty
%     errStr: If an error occured, errStr is identify the error
%
%
% Bjoern W. Kreher
% 08/04
%
% UNIX

verStr= 'V1.1';

res= [];    errStr= '';

dtdStruct= [];

if length(varargin) < 2
    errStr= 'dtdstruct_init: There have to be at least one parameter';
    return;
end

if (~isempty(varargin)) && ischar(varargin{1})
    typeStr= varargin{1};
else
    errStr= 'dtdstruct_init: First param have to be a string';
    return;
end

argNo= nargin - 1;
%argNo= argNo - 1;
param= cell(100, 1);
for i= 1:100
    if length(varargin) < (i + 1)
        param{i}= [];
    else
        param{i}= varargin{i + 1};
    end
end


if strcmp(typeStr, 'DTD')
    %Eigenvector
    if ~mrstruct_istype(param{1}) || (mrstruct_query(param{1}, 'dimensions') ~= 5)
        errStr= strcat(mfilename, ' DTD(error): Second parameter have to be a 5D mrstruct for type ');
        return
    end
    sizeEigVec= mrstruct_query(param{1}, 'sizeAy');
    if ~isequal(sizeEigVec(4:5), [3 3])
        errStr= strcat(mfilename, ' DTD(error): Second parameter have to be a 5D mrstruct witch XxYxZx3x3');
        return
    end
    dtdStruct.eigenVec_struc= param{1};
    
    %eigenvalues
    if ~mrstruct_istype(param{2}) || (mrstruct_query(param{2}, 'dimensions') ~= 4)
        errStr= strcat(mfilename, ' DTD(error): Third parameter have to be a 4D mrstruct');
        return
    end
    sizeEigVal= mrstruct_query(param{2}, 'sizeAy');
    if ~isequal(sizeEigVal(1:3), sizeEigVec(1:3)) || (sizeEigVal(4) ~= 3)
        errStr= strcat(mfilename, ' DTD(error): Third parameter have to be a 4D mrstruct witch XxYxZx3');
        return
    end
    dtdStruct.eigenVal_struc= param{2};
    dtdStruct.version= verStr;

    [res, errStr]= local_append_data(dtdStruct, param(3:2:argNo), param(4:2:argNo));
elseif strcmp(typeStr, 'VelVD')
   
    %eigenvalues
    if ~mrstruct_istype(param{1}) || (mrstruct_query(param{2}, 'dimensions') ~= 4)
        errStr= strcat(mfilename, ' VelVD(error): Third parameter have to be a 4D mrstruct');
        return
    end
    sizeEigVal= mrstruct_query(param{1}, 'sizeAy');
    if (sizeEigVal(4) ~= 3)
        errStr= strcat(mfilename, ' VelVD(error): Third parameter have to be a 4D mrstruct witch XxYxZx3');
        return
    end
    dtdStruct.velocityVect_struc= param{1};
    dtdStruct.version= verStr;

    [res, errStr]= local_append_data(dtdStruct, param(2:2:argNo), param(3:2:argNo));

elseif strcmp(typeStr, 'MR')   
    %eigenvalues
    if ~mrstruct_istype(param{1}) || (mrstruct_query(param{1}, 'dimensions') ~= 3)
        errStr= strcat(mfilename, ' MR(error): Second parameter have to be a 3D mrstruct');
        return
    end
    dtdStruct.version= verStr;
    if isempty(param{2})
        dtdStruct.mrStruct_001= param{1};
    else
        try
            dtdStruct.(param{2}) = param{1};
        catch
           errStr=  strcat(mfilename, ' MR(error): Third parameter is an invalid string A');
           return;
        end
        if strcmp(param{2}, 'version')
           errStr=  strcat(mfilename, ' MR(error): Third parameter is an invalid string');
           return;
        end            
    end

    [res, errStr]= local_append_data(dtdStruct, param(3:2:argNo), param(4:2:argNo));
    
else
    errStr= sprintf('dtdstruct_init: dtdStruct type ''%s'' is not supportrd yet', typeStr);
    return;
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_append_data(dtd, dataMR, nameCell)
res= []; errStr= '';

if length(dataMR) ~= length(nameCell)
    errStr= strcat(mfilename, '::local_append_data DTD(error): not enough names');
    return
end
for i= 1:length(dataMR)
    [dtd, errStr]= dtdstruct_modify(dtd, 'addStruct', dataMR{i}, nameCell{i}, 1);
    if ~isempty(errStr)
        return
    end
end
res= dtd;
return


