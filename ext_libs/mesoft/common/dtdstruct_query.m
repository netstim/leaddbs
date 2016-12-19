function [res, errStr, oArg]= dtdstruct_query(varargin)
%function [res, errStr]= dtdstruct_query(dtdStruct, command[, op1[, op2[... [, opN]]]])
%
%
%   General method query a infomations and data dtdStruct
%
%   dtdStruct: should contain a valid dtdStruct
%   command:    {'sizeAy' | 'vox' | 'user' | 'hMatrix' | 'patient' | 'dtdType' | 'existField' | 'getComStr' | 'getVer' |
%                'getFA' | 'getVR' | 'getRA' | 'getLinearM' | 'getPlanarM' | 'getTrace' | 'getEigVal1' | 'getEigVal2' | 'getEigVal3' |
%                'getEigVal' | 'getEigVec' | 'getDiffDir' | 'get1Dir' | 'get2Dir' | 'get3Dir' | 'getColor'}
%   op[1..end] the parameter specified by the command string
%
% return values:
%   return: the modified dtdStruct or if an error occured the old one
%   errStr: a message string, specifying the error which occured
%   oArg: additional output argument, depending on the command
% 
% Description of the different commands: 
% 'sizeAy'    returns the size of the dtdStruct
%             res: [1x3] double vector containing the size of the dtdStruct [sizeX sizeY sizeZ]
% 'vox'       returns the vox size of the dtdStruct
%             res: [1x4] double vector containing the voxel size of the dtdStruct [voxX voxY voxZ gapZ]
% 'user'      returns a the userentry identified by the sring op1
%             op1: string which identifies the user entry
%             res: user entry
% 'hMatrix'   returns the homogene transformation to transfor voxel inices to the world coord
%             res: [4x4] double matrix
% 'patient'   returns the name r identifyer if the pationt or subject
%             res: string containing the patient name
% 'dtdType'   returns the identifyer of the dtdStruct type
%             res: string to identify the type of the dtdStruct
% 'existField'tests if a given map entry name already exists
%             res: 1 if the name already exists, 0 else
% 'getComStr' determine the identifier of all available maps
%             res: cell string containing for each available map a entry string with the identifier
% 'getVer'    determines the version code of the given dtdstruct
%             res: string to identify the dtdstruct version 
% 'getFA'     determines the fractional anisotropy map
%             res: mrStruct of size [sizeX sizeY sizeZ] containing the requested map
% 'getVR'     determines the volume ratio map
%             res: mrStruct of size [sizeX sizeY sizeZ] containing the requested map
% 'getRA'     determines the relative anisotropy map
%             res: mrStruct of size [sizeX sizeY sizeZ] containing the requested map
% 'getLinearM'determines the linear anisotropy map  
%             res: mrStruct of size [sizeX sizeY sizeZ] containing the requested map
% 'getPlanarM'determines the planare anisotropy map
%             res: mrStruct of size [sizeX sizeY sizeZ] containing the requested map
% 'getTrace'  determines the mean diffusivity map
%             res: mrStruct of size [sizeX sizeY sizeZ] containing the requested map
% 'getEigVal1'determines the diffusion map of the highes eigenvalue
%             res: mrStruct of size [sizeX sizeY sizeZ] containing the requested map
% 'getEigVal2'determines the diffusion map of the middle eigenvalue  
%             res: mrStruct of size [sizeX sizeY sizeZ] containing the requested map
% 'getEigVal3'determines the diffusion map of the lowest eigenvalue
%             res: mrStruct of size [sizeX sizeY sizeZ] containing the requested map
% 'getXXX'    determines the map specified by the the identifier (by replacing XXX by the identifyer).
%             the identifyer can be determined by the 'getComStr' command
%             res: mrStruct of size [sizeX sizeY sizeZ] containing the requested map            
% 'getEigVal' returns all three eigenvalues
%             res: mrStruct of size [sizeX sizeY sizeZ 3] containing the requested map
% 'getEigVec' returns all three eigenvectors 
%             res: mrStruct of size [sizeX sizeY sizeZ 3 3] containing the requested map
% 'getDiffDir'returns the direction of highest difusivity (equivalent ti get1Dir)
%             res: mrStruct of size [sizeX sizeY sizeZ 3] containing the requested vector field
% 'get1Dir'   returns the eigenvector corresponding to the highest eigenvalue
%             res: mrStruct of size [sizeX sizeY sizeZ 3] containing the requested vector field
% 'get2Dir'   returns the eigenvector corresponding to the middle eigenvalue
%             res: mrStruct of size [sizeX sizeY sizeZ 3] containing the requested vector field
% 'get3Dir'   returns the eigenvector corresponding to the lowest eigenvalue
%             res: mrStruct of size [sizeX sizeY sizeZ 3] containing the requested vector field
% 'getColor'  determines the unscaled color coded directionmap
%             res: int8 array of the size [sizeX sizeY sizeZ 3] containing the requested color map
% 
%
% Bjoern W. Kreher
% 07/02
%
% UNIX







res= [];    errStr= '';     oArg= [];

[isT, dtdType, errStr]= dtdstruct_istype(varargin{1});
if isT && strcmp(dtdType, 'mrStruct')
    varargin{1}= dtdstruct_init('MR', varargin{1});
end

if ischar(varargin{1}) && exist(varargin{1}, 2)
    dtdStruct= dtdstruct_read(varargin{1});
elseif isT
    dtdStruct= varargin{1};
else
    errStr= 'Error in dtdstruct_query: first param have to be dtdStruct';
    return;
end

if ~ischar(varargin{2})
    errStr= 'Error in dtdstruct_query: command have to be a string';
    return;
end

for i= 1:7
    if length(varargin) < (i + 2)
        param{i}= [];
    else
        param{i}= varargin{i + 2};
    end
end

commandStr= varargin{2};

version= local_getVer(dtdStruct);

mrStruComStr= {'_mrStruct'};
%dtdComStr= {'FA'; 'Trace'; 'EigVal1'; 'EigVal2'; 'EigVal3'};
dtdComStr= {'FA'; 'Trace'; 'EigVal1'; 'EigVal2'; 'EigVal3'; 'VR'; 'RA'; 'LinearM'; 'PlanarM'; 'SphericalM'};
monoMdtdComStr= {'FA_1'; 'TrD_1'; 'TrD_ISO'; 'Fraction1'; 'FractionISO'};
mdtdComStr= {'FA_2'; 'TrD_2'; 'Fraction2'};
vdtComStr= {'Magnitude'};

dtdVecComStr= {'1Dir'; '2Dir'; '3Dir'};
mdtVecComStr= {'1DirA'; '2DirA'; '3DirA'; '1DirB'; '2DirB'; '3DirB'};
monoMdtVecComStr= {'1DirA'; '2DirA'; '3DirA'};
vdtVecComStr= {'Dir'};

if strcmp(commandStr, 'sizeAy')
    [res, errStr]= local_getSizeAy(dtdStruct);
    return;
elseif strcmp(commandStr, 'vox')
    dtdCell= struct2cell(dtdStruct);
    i= 1;
    while ~mrstruct_istype(dtdCell{i})
        i= i + 1;
    end
    [res, errStr]= local_getVox(dtdCell{i});
    return;
elseif strcmp(commandStr, 'user')
    if isempty(param{1}) || ischar(param{1})
        userStr= param{1};
    else
        errStr= strcat(mfilename, '(error): Paramter userStr should be empty or of type string');
        return 
    end
    [res, errStr]= local_getUser(dtdStruct, userStr);
    return;
elseif strcmp(commandStr, 'hMatrix')
    [res, errStr]= local_getHMatrix(dtdStruct);
    return;
elseif strcmp(commandStr, 'patient')
    dtdCell= struct2cell(dtdStruct);
    i= 1;
    while ~mrstruct_istype(dtdCell{i})
        i= i + 1;
    end
    res= dtdCell{i}.patient;
    return;
    
elseif strcmp(commandStr, 'cmpDimensions')
    dtdCell= struct2cell(dtdStruct);
    for i= 1:length(dtdCell)
        if mrstruct_istype(dtdCell{i})
            [res, errStr]= mrstruct_query(dtdCell{i}, 'cmpDimensions', param{1});
            if ~isempty(res) && (res > 0)
                res= true;  errStr= '';
            end
            return;
        end
    end
    errStr= strcat(mfilename, '(error): Internal error (struct contains no mr)');
    return;
elseif strcmp(commandStr, 'mrStructProb')
    dtdCell= struct2cell(dtdStruct);
    for i= 1:length(dtdCell)
        if mrstruct_istype(dtdCell{i})
            res= mrstruct_init('volume', [], dtdCell{i});
            return;
        end
    end
    errStr= strcat(mfilename, '(error): Internal error (struct contains no mr)');
    return;
elseif strcmp(commandStr, 'dtdType')
    res= dtdType;
    return;
elseif strcmp(commandStr, 'existField')
    [res, errStr]= local_existField(dtdStruct, param{1});
    return;
elseif strcmp(commandStr, 'getComStr')
    if strcmp(dtdType, 'mrStruct')
        res= mrStruComStr;
    else
        [res, errStr]= local_getComStr(dtdStruct, dtdType, dtdComStr, mdtdComStr, monoMdtdComStr, vdtComStr, 3);
        [oArg, errStr]= local_getComStr(dtdStruct, dtdType, dtdVecComStr, mdtVecComStr, monoMdtVecComStr, vdtVecComStr, 4);
    end
    return;
elseif strcmp(commandStr, 'getVer')
    res= version;   errStr= '';
    return;
end
if strcmp(dtdType, 'VelVD')
    if strcmp(commandStr, 'getMagnitude')
        [res, errStr]= local_getMagnitude(dtdStruct, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, 'getDiffDir') || strcmp(commandStr, 'getDir')
        res= local_getVelocityDir(dtdStruct.velocityVect_struc, 'no', param{1}, param{2}, param{3});
        res= local_coordCorr(res, version, param{4});
        return; 
    end
elseif strcmp(dtdType, 'multiDTD') || strcmp(dtdType, 'DTD') || strcmp(dtdType, 'monoMDT')
    if strcmp(commandStr, 'getFA') % getFA
        [res, errStr]= local_getFA(dtdStruct, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, 'getVR') % getTrace
        [res, errStr]= local_getVR(dtdStruct, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, 'getRA') % getTrace
        [res, errStr]= local_getRA(dtdStruct, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, 'getLinearM') % getLinearM (see Alexander et al: CL)
        [res, errStr]= local_getLinearM(dtdStruct, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, 'getPlanarM') % getPlanarM (see Alexander et al: CP)
        [res, errStr]= local_getPlanarM(dtdStruct, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, 'getSphericalM') % getSpericalM (see Alexander et al: CS)
        [res, errStr]= local_sphericM(dtdStruct, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, strcat('get', 'Trace')) % getTrace
        [res, errStr]= local_getTrace(dtdStruct, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, 'getEigVal1')
        [res, errStr]= local_getEigValSort(dtdStruct.eigenVal_struc, 1, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, 'getEigVal2')
        [res, errStr]= local_getEigValSort(dtdStruct.eigenVal_struc, 2, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, 'getEigVal3')
        [res, errStr]= local_getEigValSort(dtdStruct.eigenVal_struc, 3, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, 'getEigVal')
        [res, errStr]= local_getEigVal(dtdStruct.eigenVal_struc, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, 'getEigVec')
        [res, errStr]= local_getEigVec(dtdStruct.eigenVec_struc, param{1}, param{2}, param{3});
        res= local_coordCorr(res, version, '');
        return;
        
    elseif strcmp(commandStr, 'getDiffDir') || strcmp(commandStr, 'get1Dir')
        [res, errStr]= local_getDiffDir(dtdStruct.eigenVal_struc, dtdStruct.eigenVec_struc, 1, param{1}, param{2}, param{3});
        res= local_coordCorr(res, version, param{4});
        return;
    elseif strcmp(commandStr, 'get2Dir')
        [res, errStr]= local_getDiffDir(dtdStruct.eigenVal_struc, dtdStruct.eigenVec_struc, 2, param{1}, param{2}, param{3});
        res= local_coordCorr(res, version, param{4});
        return;
    elseif strcmp(commandStr, 'get3Dir')
        [res, errStr]= local_getDiffDir(dtdStruct.eigenVal_struc, dtdStruct.eigenVec_struc, 3, param{1}, param{2}, param{3});
        res= local_coordCorr(res, version, param{4});
        return;
    elseif strcmp(commandStr, 'getColor')
        [res, errStr]= local_getColor(dtdStruct.eigenVal_struc, dtdStruct.eigenVec_struc, param{1}, param{2}, param{3}, param{4}, param{5}, version);    
        %    res= res(:, :, :, [2 1 3]);
        return;
    elseif strcmp(commandStr, 'calcDW')
        [res, errStr]= local_getDWImages(dtdStruct.eigenVal_struc, dtdStruct.eigenVec_struc, param{1}, version);            
    end
end
if strcmp(dtdType, 'mrStruct')
    
end
if strcmp(dtdType, 'multiDTD') || strcmp(dtdType, 'monoMDT')
    if strcmp(commandStr, 'getDiffDirA') || strcmp(commandStr, 'get1DirA')
        [res, errStr]= local_getDiffDir(dtdStruct.eigVal_1, dtdStruct.eigVec_1, 1, param{1}, param{2}, param{3});
        res= local_coordCorr(res, version, param{4});
        return;
    elseif strcmp(commandStr, 'getEigVal_A')
        [res, errStr]= local_getEigVal(dtdStruct.eigVal_1, param{1}, param{2}, param{3});
    elseif strcmp(commandStr, 'getEigVec_A')
        [res, errStr]= local_getEigVec(dtdStruct.eigVec_1, param{1}, param{2}, param{3});
        res= local_coordCorr(res, version, '');
    elseif strcmp(commandStr, 'get2DirA')
        [res, errStr]= local_getDiffDir(dtdStruct.eigVal_1, dtdStruct.eigVec_1, 2, param{1}, param{2}, param{3});
        res= local_coordCorr(res, version, param{4});
        return;
    elseif strcmp(commandStr, 'get3DirA')
        [res, errStr]= local_getDiffDir(dtdStruct.eigVal_1, dtdStruct.eigVec_1, 3, param{1}, param{2}, param{3});
        res= local_coordCorr(res, version, param{4});
        return;
    elseif strcmp(commandStr, 'getColorA')
        [res, errStr]= local_getColor(dtdStruct.eigVal_1, dtdStruct.eigVec_1, param{1}, param{2}, param{3}, param{4}, param{5}, version);
        return;
        
    elseif strcmp(commandStr, strcat('get', monoMdtdComStr{1})) % getFA_1
        [res, errStr]= local_calcFA(dtdStruct.eigVal_1, [], param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, strcat('get', monoMdtdComStr{2})) % getTrD_1
        [res, errStr]= local_calcTrD(dtdStruct.eigVal_1, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, strcat('get', monoMdtdComStr{3})) % getTrD_ISO
        [res, errStr]= local_calcTrD(dtdStruct.eigVal_3, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, 'getFraction1')
        [res, errStr]= local_getFraction(dtdStruct, 1, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, 'getFractionISO')
        if strcmp(dtdType, 'monoMDT')
            [res, errStr]= local_getFraction(dtdStruct, 2, param{1}, param{2}, param{3});
        else
            [res, errStr]= local_getFraction(dtdStruct, 3, param{1}, param{2}, param{3});
        end            
        return;
    end
end
if strcmp(dtdType, 'multiDTD')  
    if strcmp(commandStr, 'getDiffDirB') || strcmp(commandStr, 'get1DirB')
        [res, errStr]= local_getDiffDir(dtdStruct.eigVal_2, dtdStruct.eigVec_2, 1, param{1}, param{2}, param{3});
        res= local_coordCorr(res, version, param{4});
        return;
    elseif strcmp(commandStr, 'getEigVal_B')
        [res, errStr]= local_getEigVal(dtdStruct.eigVal_2, param{1}, param{2}, param{3});
    elseif strcmp(commandStr, 'getEigVec_B')
        [res, errStr]= local_getEigVec(dtdStruct.eigVec_2, param{1}, param{2}, param{3});
        res= local_coordCorr(res, version, '');
    elseif strcmp(commandStr, 'get2DirB')
        [res, errStr]= local_getDiffDir(dtdStruct.eigVal_2, dtdStruct.eigVec_2, 2, param{1}, param{2}, param{3});
        res= local_coordCorr(res, version, param{4});
        return;
    elseif strcmp(commandStr, 'get3DirB')
        [res, errStr]= local_getDiffDir(dtdStruct.eigVal_2, dtdStruct.eigVec_2, 3, param{1}, param{2}, param{3});
        res= local_coordCorr(res, version, param{4});
        return;
    elseif strcmp(commandStr, 'getColorB')
        [res, errStr]= local_getColor(dtdStruct.eigVal_2, dtdStruct.eigVec_2, param{1}, param{2}, param{3}, param{4}, param{5}, version);
        return;
    elseif strcmp(commandStr, strcat('get', mdtdComStr{1})) % getFA_2
        [res, errStr]= local_calcFA(dtdStruct.eigVal_2, [], param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, strcat('get', mdtdComStr{2})) % getTrD_2
        [res, errStr]= local_calcTrD(dtdStruct.eigVal_2, param{1}, param{2}, param{3});
        return;
    elseif strcmp(commandStr, 'getFraction2')
        [res, errStr]= local_getFraction(dtdStruct, 2, param{1}, param{2}, param{3});
        return;
    end
end


if (length(commandStr) >= 4) && strcmp(commandStr(1:4), 'get_') % generated Data
    [res, errStr]= local_getGeneral(dtdStruct, commandStr(5:end), param{1}, param{2}, param{3});
    return;
end

errStr= sprintf('Error in dtdstruct_query: command ''%s'' is not implemented for dtd type ''%s''', commandStr, dtdType);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%
%  START:
%       [result, errStr]= local_getVer(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result, errStr]= local_getVer(dtdStruct)

errStr= '';
if isfield(dtdStruct, 'version')
    result= dtdStruct.version;
else
    result= 'V1.0'; 
end


%
%
%  START:
%       function [res, errStr]= local_getUser(dtdStruct, userStr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_getUser(dtdStruct, userStr)
res= []; errStr= '';

dtdCell= struct2cell(dtdStruct);
i= 1;
while ~mrstruct_istype(dtdCell{i})
    i= i + 1;
    if i < length(dtdCell)
        errStr= strcat(mfilename, ' (error): dtdStruct contains mrStruct');
        return
    end
end
userStrc= dtdCell{i}.user;

if isempty(userStr)
    res= userStrc;
else
    if ~isempty(userStrc)
       userCell= struct2cell(userStrc);
       userNames= fieldnames(userStrc);
       idx= find(strcmp(userNames, userStr));
    else
        idx= [];
    end
   if isempty(idx)
       errStr= strcat(mfilename, '::local_getUser(error): field name ''', userStr, ''' was not found in dtdStruct');
       return
   end
   res= userCell{idx};
end

%
%
%  START:
%       function [res, errStr]= local_getSizeAy(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getSizeAy(dtdStruct)
res= []; errStr= '';

dtdCell= struct2cell(dtdStruct);
i= 1;
while ~mrstruct_istype(dtdCell{i})
    i= i + 1;
    if i > length(dtdCell)
        errStr= strcat(mfilename, '::local_getSizeAy(error): dtdStruct contains no valid mrStruct');
        return
    end
end
sizeAy= mrstruct_query(dtdCell{i}, 'sizeAy');
if length(sizeAy) < 3
    res= ones(1, 3);
    res(1:length(sizeAy))= sizeAy;
else
    res= sizeAy(1:3);
end

%
%
%  START:
%       function [res, errStr]= local_getHMatrix(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getHMatrix(dtdStruct)
res= []; errStr= '';

dtdCell= struct2cell(dtdStruct);
i= 1;
while ~mrstruct_istype(dtdCell{i})
    i= i + 1;
    if i > length(dtdCell)
        errStr= strcat(mfilename, '::local_getHMatrix(error): dtdStruct contains no valid mrStruct');
        return
    end
end

if isnumeric(dtdCell{i}.edges) && isequal(size(dtdCell{i}.edges), [4 4])
    res= dtdCell{i}.edges;
elseif isfield(dtdCell{i}.user, 'hMatrix')
    res= dtdCell{i}.user.hMatrix;
else
    errStr= strcat(mfilename, '::local_getHMatrix(warning): dtdStruct contains no valid transformation matrix');    
end

%
%
%  START:
%       [result, errStr]= local_coordCorr(mrDir, version, coord_flag)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mrDir, errStr]= local_coordCorr(mrDir, version, coord_flag)

if strcmp(version, 'V1.0') && isempty(coord_flag)
    if length(size(mrDir.dataAy)) == 4 % vectors
        mrDir.dataAy= mrDir.dataAy(:, :, :, [2 1 3]);
        mrDir.dataAy(:, :, :, 3)= -mrDir.dataAy(:, :, :, 3);
    elseif length(size(mrDir.dataAy)) == 5  % rotMatrix
        mrDir.dataAy= mrDir.dataAy(:, :, :, [2 1 3], :);
        mrDir.dataAy(:, :, :, 3, :)= -mrDir.dataAy(:, :, :, 3, :);
    end
end

        
%
%
%  START:
%       [result, errStr]= local_getEigenVec(eigenVec, x, y, z)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getDWImages(eigVal, eigVec, bTen, version)
res= []; errStr= '';

sizeAy= mrstruct_query(eigVec, 'sizeAy');
bTenNo= size(bTen, 3);
res= mrstruct_init('series3D', zeros([sizeAy(1:3) bTenNo]), eigVal);

mask= sum(abs(eigVal.dataAy), 4) > 0;
bTenR= reshape(bTen, [9 size(bTen, 3)]);
for z= 1:sizeAy(3)
    for x= 1:sizeAy(1)
        for y= 1:sizeAy(2)
            if mask(x, y, z) ~= 0
                eVec= squeeze(eigVec.dataAy(x, y, z, :, :));
                eVal= diag(squeeze(eigVal.dataAy(x, y, z, :)));
                dti= reshape(eVec*eVal*(eVec'), [9 1])*ones(1, size(bTenR, 2));                
                res.dataAy(x, y, z, :)= exp(-sum(dti.*bTenR, 1));                
            end
        end
    end
    disp(z);
end

res.user.B_Tensor= bTen;

%
%
%  START:
%       [result, errStr]= local_getEigenVec(eigenVec, x, y, z)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getEigVec(eigenVec, x, y, z)
errStr= '';
[x, y, z]= local_getArea(size(eigenVec.dataAy), x, y, z);

res= mrstruct_init('series3DEchos', eigenVec.dataAy(x, y, z, :, :), eigenVec);

%res= eigenVec;
%res.dataAy= eigenVec.dataAy(x, y, z, :, :);

%
%
%  START:
%       [result, errStr]= local_getEigenValSort(eigenVec, no, x, y, z)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getEigValSort(eigenVal, no, x, y, z)
res= []; errStr= '';

if isempty(no) || ~isnumeric(no) || (numel(no) ~= 1) || ((1 ~= no) && (2 ~= no) && (3 ~= no))
    errStr= strcat(mfilename, ':local_getEigValSort (error): no have to bi an single integer between 1 and three');
    return
end
[x, y, z]= local_getArea(size(eigenVal.dataAy), x, y, z);

tmpAy= sort(eigenVal.dataAy(x, y, z, :), 4);
res= mrstruct_init('volume', tmpAy(:, :, :, 4 - no), eigenVal);

%
%
%  START:
%       [result, errStr]= local_getEigenVal(eigenVec, x, y, z)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getEigVal(eigenVal, x, y, z)
errStr= '';
[x, y, z]= local_getArea(size(eigenVal.dataAy), x, y, z);

res= mrstruct_init('series3D', eigenVal.dataAy(x, y, z, :), eigenVal);


%
%
%  START:
%       [result, errStr]= local_getVelocityDir(velVec, normFlag, x, y, z)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result, errStr]= local_getVelocityDir(velVec, normFlag, x, y, z)

errStr= '';
[x, y, z]= local_getArea(size(velVec.dataAy), x, y, z);

if isempty(normFlag) || ~strcmp(normFlag, 'yes')
    result= mrstruct_init('series3D', velVec.dataAy(x, y, z, :), velVec);    
else
    tmpMy= velVec.dataAy(x, y, z, :);
    lenMy= sum(tmpMy.^2, 4);
    idx= find(lenMy > 0);
    lenMy(idx)= 1./sqrt(lenMy(idx));
    tmpMy(:, :, :, 1)= tmpMy(:, :, :, 1).*lenMy;
    tmpMy(:, :, :, 2)= tmpMy(:, :, :, 2).*lenMy;
    tmpMy(:, :, :, 3)= tmpMy(:, :, :, 3).*lenMy;
    result= mrstruct_init('series3D', tmpMy, velVec);
end

%
%
%  START:
%       [result, errStr]= local_getDiffDir(eigenVal, eigenVec)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result, errStr]= local_getDiffDir(eigenVal, eigenVec, eigNo, x, y, z)

errStr= '';
[x, y, z]= local_getArea(size(eigenVal.dataAy), x, y, z);

eigenVal.dataAy= eigenVal.dataAy(x, y, z, :);
eigenVec.dataAy= eigenVec.dataAy(x, y, z, :, :);

sVal=size(eigenVal.dataAy);
sVal= sVal(1:4);
result= mrstruct_init('series3D', [], eigenVal);


if isstruct(eigenVec.user) && (sum(strcmp(fieldnames(eigenVec.user), 'sort')) > 0)
    order= eigenVec.user.sort;
    result.dataAy= reshape(eigenVec.dataAy(:, :, :, :, order(eigNo)), sVal);
else
    %perEigVec= permute(squeeze(eigenVec.dataAy), [5 1 2 3 4]);
    perEigVec= permute((eigenVec.dataAy), [4 1 2 3 5]);
    temp= reshape(eigenVal.dataAy, [numel(eigenVal.dataAy)/3 3]);
    [dummy, index]= sort(abs(temp'));
    [dummy, index2]= max(abs(temp'));
    clear dummy;
    clear temp;
    ss= (index(4 - eigNo, :) - 1)*numel(perEigVec)/3 + (1:3:numel(perEigVec)/3);
    temp= reshape(perEigVec, [numel(perEigVec)/3 3]);
    
    sss= [ss ss+1 ss+2]; % <- entspricht permute
    result.dataAy= reshape(temp(sss), sVal);
end

%
%
%  START:
%       function [res, errStr]= local_getColor(dirVc, backGround, mask)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [imageRGB, errStr]= local_getColor(eigVal, eigVec, x, y, z, imageRGB, mask, ver)
errStr= '';

[x, y, z]= local_getArea(size(eigVal.dataAy), x, y, z);

[dirAy, errStr]= local_getDiffDir(eigVal, eigVec, 1, x, y, z);
dirAy= local_coordCorr(dirAy, ver, '');

if isempty(imageRGB) || isempty(mask)
    imageRGB= uint8(255*(1 - 2*acos(abs(dirAy.dataAy))/pi));
    imageRGB= imageRGB(:, :, :, [2 1 3]);
    return;
end

base= numel(dirAy.dataAy)/3;
index= find(mask ~= 0);
perm_comp= [2 1 3];
dirRGB= uint8(zeros(length(index), 3));
alphaRGB= [1 1 1];
alphaScale= 1;
transRGB= [[1   0   0]; ...
           [0   0.8 0]; ...
           [0   0   1]];
for i=  1:3
    if 1 %skalierter RGB-Wert mit alpha
        scale= 255*((double(imageRGB((i - 1)*base + index))/255).^(alphaRGB(i)*alphaScale));
        dirRGB(:, i)= uint8(floor(scale.*(1 - 2*acos(abs(dirAy.dataAy((perm_comp(i) - 1)*base + index)))/pi)));
    elseif 1 % skalierter RGB-Wert mit farb transformation (transRGB)
        tmpRGB= zeros(length(index), 1);
        for ii= 1:3
            tmpRGB(:)= tmpRGB(:) + transRGB(i, ii)*(1 - 2*acos(abs(dirAy.dataAy((perm_comp(ii) - 1)*base + index)))/pi);
        end
        dirRGB(:, i)= uint8(floor((255*255^-alphaScale)*(double(imageRGB((i - 1)*base + index)).^alphaScale).*tmpRGB));
    else % unskaliertrer RGB wert
        dirRGB(:, i)= uint8(floor(255*(1 - 2*acos(abs(dirAy.dataAy((perm_comp(i) - 1)*base + index)))/pi)));
    end
end
%dirRGB= dirRGB(:, [2 1 3]);



imageRGB(0*base + index)= uint8(dirRGB(:, 1));
imageRGB(1*base + index)= uint8(dirRGB(:, 2));
imageRGB(2*base + index)= uint8(dirRGB(:, 3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getColor
%

%
%
%  START:
%       function [res, errStr]= local_existField(dtdStruct, fieldNameStr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xAy, yAy, zAy]= local_getArea(szAy, x, y, z)

sizeAy= ones(1, 4);
sizeAy(1:length(szAy))= szAy;
if isempty(x),  xAy= 1:sizeAy(1);
else            xAy= x;
end
if isempty(y),  yAy= 1:sizeAy(2);
else            yAy= y;
end
if isempty(z),  zAy= 1:sizeAy(3);
else            zAy= z;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_existField
%


%
%
%  START:
%       function [res, errStr]= local_existField(dtdStruct, fieldNameStr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_existField(dtdStruct, fieldNameStr)

errStr= '';
res= 0;

stNames= fieldnames(dtdStruct);
res= sum(strcmp(stNames, fieldNameStr));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_existField
%

%
%
%  START:
%       [res, errStr]= local_getComStr(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getComStr(dtdStruct, typeStr, dtdComStr, mdtdComStr, monoMDTComStr, vdtComStr, dimNo)

errStr= '';
res= {};
stNames= fieldnames(dtdStruct);
dtdCell= struct2cell(dtdStruct);
if strcmp(typeStr, 'DTD')
    res= dtdComStr;
elseif strcmp(typeStr, 'monoMDT')
    res= {dtdComStr{:} monoMDTComStr{:}};
elseif strcmp(typeStr, 'multiDTD')
    res= {dtdComStr{:} monoMDTComStr{:} mdtdComStr{:}};
elseif strcmp(typeStr, 'VelVD')
    res= vdtComStr;
end

for i= 1:length(dtdCell)
    if mrstruct_istype(dtdCell{i}) && length(size(dtdCell{i}.dataAy)) == dimNo
        if ~strcmp(stNames{i}, 'TrD_image_struc') && ~strcmp(stNames{i}, 'FA_image_struc') ...
                && ~strcmp(stNames{i}, 'eigenVal_struc') && ~strcmp(stNames{i}, 'eigVal_1') && ~strcmp(stNames{i}, 'eigVal_2') ...
                && ~strcmp(stNames{i}, 'velocityVect_struc')
            res{end+1, 1}= strcat('_', stNames{i});
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getTrace
%
%

%
%  START:
%       [res, errStr]= local_getLinearM(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getLinearM(dtdStruct, x, y, z)

res= []; errStr= '';
[x, y, z]= local_getArea(local_getSizeAy(dtdStruct), x, y, z);
sizeAy= [length(x), length(y), length(z)];

[eVal1, errStr]= local_getEigValSort(dtdStruct.eigenVal_struc, 1, x, y, z);
[eVal2, errStr]= local_getEigValSort(dtdStruct.eigenVal_struc, 2, x, y, z);
[eVal3, errStr]= local_getEigValSort(dtdStruct.eigenVal_struc, 3, x, y, z);

inDataAy= dtdStruct.eigenVal_struc.dataAy(x, y, z, :);
idx= find(min(inDataAy, [], 4) > 0);
dataAy= zeros(sizeAy);
sumAy= (inDataAy(0*prod(sizeAy) + idx) + inDataAy(1*prod(sizeAy) + idx) + inDataAy(2*prod(sizeAy) + idx));

dataAy(idx)= (eVal1.dataAy(idx) - eVal2.dataAy(idx))./sumAy;

res= mrstruct_init('volume', dataAy, dtdStruct.eigenVal_struc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getLinearM
%

%
%  START:
%       [res, errStr]= local_getLinearM(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getPlanarM(dtdStruct, x, y, z)

res= []; errStr= '';
[x, y, z]= local_getArea(local_getSizeAy(dtdStruct), x, y, z);
sizeAy= [length(x), length(y), length(z)];

[eVal1, errStr]= local_getEigValSort(dtdStruct.eigenVal_struc, 1, x, y, z);
[eVal2, errStr]= local_getEigValSort(dtdStruct.eigenVal_struc, 2, x, y, z);
[eVal3, errStr]= local_getEigValSort(dtdStruct.eigenVal_struc, 3, x, y, z);

inDataAy= dtdStruct.eigenVal_struc.dataAy(x, y, z, :);
idx= find(min(inDataAy, [], 4) > 0);
dataAy= zeros(sizeAy);
sumAy= (inDataAy(0*prod(sizeAy) + idx) + inDataAy(1*prod(sizeAy) + idx) + inDataAy(2*prod(sizeAy) + idx));

dataAy(idx)= 2*(eVal2.dataAy(idx) - eVal3.dataAy(idx))./sumAy;

res= mrstruct_init('volume', dataAy, dtdStruct.eigenVal_struc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getLinearM
%


%
%  START:
%       [res, errStr]= local_getRA(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getRA(dtdStruct, x, y, z)

res= []; errStr= '';
[x, y, z]= local_getArea(local_getSizeAy(dtdStruct), x, y, z);
sizeAy= [length(x), length(y), length(z)];
inDataAy= dtdStruct.eigenVal_struc.dataAy(x, y, z, :);
idx= find(min(inDataAy, [], 4) > 0);
dataAy= zeros(sizeAy);
meanAy= (inDataAy(0*prod(sizeAy) + idx) + inDataAy(1*prod(sizeAy) + idx) + inDataAy(2*prod(sizeAy) + idx))./3;
tmp1Ay= (inDataAy(0*prod(sizeAy) + idx) - meanAy).^2 + (inDataAy(1*prod(sizeAy) + idx) - meanAy).^2 + (inDataAy(2*prod(sizeAy) + idx) - meanAy).^2;
dataAy(idx)= sqrt(tmp1Ay./(3*meanAy));
res= mrstruct_init('volume', dataAy, dtdStruct.eigenVal_struc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getRA
%

%
%  START:
%       [res, errStr]= local_getVR(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getVR(dtdStruct, x, y, z)

res= []; errStr= '';
[x, y, z]= local_getArea(local_getSizeAy(dtdStruct), x, y, z);
sizeAy= [length(x), length(y), length(z)];
inDataAy= dtdStruct.eigenVal_struc.dataAy(x, y, z, :);
idx= find(min(inDataAy, [], 4) > 0);
dataAy= zeros(sizeAy);
tmp1Ay= inDataAy(0*prod(sizeAy) + idx) .* inDataAy(1*prod(sizeAy) + idx) .* inDataAy(2*prod(sizeAy) + idx);
meanAy= (inDataAy(0*prod(sizeAy) + idx) + inDataAy(1*prod(sizeAy) + idx) + inDataAy(2*prod(sizeAy) + idx))./3;
dataAy(idx)= tmp1Ay./(meanAy.^3);
res= mrstruct_init('volume', dataAy, dtdStruct.eigenVal_struc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getVR
%

%
%  START:
%       [res, errStr]= local_getVR(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_sphericM(dtdStruct, x, y, z)

res= []; errStr= '';
[x, y, z]= local_getArea(local_getSizeAy(dtdStruct), x, y, z);
sizeAy= [length(x), length(y), length(z)];

[eVal1, errStr]= local_getEigValSort(dtdStruct.eigenVal_struc, 1, x, y, z);
[eVal2, errStr]= local_getEigValSort(dtdStruct.eigenVal_struc, 2, x, y, z);
[eVal3, errStr]= local_getEigValSort(dtdStruct.eigenVal_struc, 3, x, y, z);

inDataAy= dtdStruct.eigenVal_struc.dataAy(x, y, z, :);
idx= find(min(inDataAy, [], 4) > 0);
dataAy= zeros(sizeAy);
mean_eVal = (eVal1.dataAy(idx) + eVal2.dataAy(idx) + eVal3.dataAy(idx))./3;
dataAy(idx)=  eVal3.dataAy(idx)./mean_eVal;
%dataAy(idx)= 3./(eVal2.dataAy(idx) + eVal3.dataAy(idx));

res= mrstruct_init('volume', dataAy, dtdStruct.eigenVal_struc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getVR
%

%
%  START:
%       [res, errStr]= local_getTrace(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getTrace(dtdStruct, x, y, z)

errStr= '';
[x, y, z]= local_getArea(local_getSizeAy(dtdStruct), x, y, z);
if dtdstruct_query(dtdStruct, 'existField', 'TrD_image_struc')
    res= dtdStruct.TrD_image_struc;
    res.dataAy= res.dataAy(x, y, z);
else
    res= local_calcTrD(dtdStruct.eigenVal_struc, x, y, z);
    res.dataAy= res.dataAy/3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getTrace
%

%
%
%  START:
%       [res, errStr]= local_getError(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getError(dtdStruct, x, y, z)

errStr= '';
[x, y, z]= local_getArea(local_getSizeAy(dtdStruct), x, y, z);

res= dtdStruct.error_struc;
res.dataAy= dtdStruct.error_struc, dataAy(x, y, z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getError
%


%
%
%  START:
%       [res, errStr]= local_getMagnitude(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getMagnitude(dtdStruct, x, y, z)

errStr= '';
[x, y, z]= local_getArea(local_getSizeAy(dtdStruct), x, y, z);

res= mrstruct_init('volume', sqrt(sum(dtdStruct.velocityVect_struc.dataAy(x, y, z, :).^2, 4)), dtdStruct.velocityVect_struc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getMagnitude
%

%
%
%  START:
%       [res, errStr]= local_getFA(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getFA(dtdStruct, x, y, z)

errStr= '';
[x, y, z]= local_getArea(local_getSizeAy(dtdStruct), x, y, z);

if dtdstruct_query(dtdStruct, 'existField', 'FA_image_struc')
    res= dtdStruct.FA_image_struc;
    res.dataAy= dtdStruct.FA_image_struc.dataAy(x, y, z);
else
    if dtdstruct_query(dtdStruct, 'existField', 'TrD_image_struc')
        [res, errStr]= local_calcFA(dtdStruct.eigenVal_struc, dtdStruct.TrD_image_struc, x, y, z);
    else
        [res, errStr]= local_calcFA(dtdStruct.eigenVal_struc, [], x, y, z);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getFA
%


%
%
%  START:
%       [res, errStr]= local_calcTrD(eigVal)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_calcTrD(eigVal, x, y, z)

errStr= '';
[x, y, z]= local_getArea(size(eigVal.dataAy), x, y, z);

res= mrstruct_init('volume', sum(eigVal.dataAy(x, y, z, :), 4), eigVal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_calcFA
%
%


%
%  START:
%       [res, errStr]= local_calcFA(eigVal, trD)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_calcFA(eigVal, trD, x, y, z)

errStr= '';
[x, y, z]= local_getArea(size(eigVal.dataAy), x, y, z);
if isempty(trD)
    trD= local_calcTrD(eigVal, x, y, z);
    sizeAy= size(trD.dataAy);
    res= trD;
else
    trD.dataAy= trD.dataAy(x, y, z);
end

eigVal.dataAy= eigVal.dataAy(x, y, z, :);
sizeAy= size(trD.dataAy);
res= trD;
res.dataAy= zeros(sizeAy);

eigV= reshape(eigVal.dataAy, [prod(sizeAy), 3]);
trD.dataAy= reshape(trD.dataAy, [prod(sizeAy), 1]);
idx= find(trD.dataAy ~= 0);
if ~isempty(idx)
    res.dataAy(idx)= 3*((eigV(idx, 1) - trD.dataAy(idx)/3).^2 + (eigV(idx, 2) - trD.dataAy(idx)/3).^2 + (eigV(idx, 3) - trD.dataAy(idx)/3).^2)/2;
    res.dataAy(idx)= sqrt(res.dataAy(idx)./sum(eigV(idx, :).^2, 2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_calcFA
%


%
%
%  START:
%       [res, errStr]= local_getVox(dtdStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getVox(mrStruct)

errStr= '';
tmp= reshape(mrStruct.vox, [1 numel(mrStruct.vox)]);
if length(tmp) == 3
    res= [tmp(1:3) 0];
elseif length(tmp) == 4
    res= tmp;
else
    res= [];
    errStr= 'dtdStruct_query(local_getVox); Warning: no valid vox size';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getVox
%


%
%
%  START:
%       [res, errStr]= local_getFraction(dtdStruct, tenId)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getFraction(dtdStruct, tenId, x, y, z)

errStr= '';
[x, y, z]= local_getArea(local_getSizeAy(dtdStruct), x, y, z);

dtdCell= struct2cell(dtdStruct);
i= 1;
while ~mrstruct_istype(dtdCell{i})
    i= i + 1;
end

res= mrstruct_init('series3D', [], dtdCell{i});

if isempty(dtdStruct.fraction)
    res.dataAy= (res.dataAy > std(reshape(res.dataAy, [numel(res.dataAy) 1]))) .* rand(size(res.dataAy)) + 10;
else
    res.dataAy= abs(dtdStruct.fraction.dataAy(x, y, z, tenId));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getFraction
%

%
%
%  START:
%       [res, errStr]= local_getGeneral(dtdStruct, strucName)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_getGeneral(dtdStruct, strucName, x, y, z)

res= []; errStr= '';

nameStr= fieldnames(dtdStruct);
dtdCell= struct2cell(dtdStruct);

idx= strcmp(strucName, nameStr);
if isempty(idx)
    errStr= strcat('dtdstruct_query::local_getGeneral (error): No mrStruct ''', strucName, ''' defined');
    return;
end
res= dtdCell{idx};
sizeAy= mrstruct_query(res, 'sizeAy');
[x, y, z]= local_getArea(sizeAy, x, y, z);
res.dataAy= real(res.dataAy(x, y, z,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  END: local_getGeneral
%

