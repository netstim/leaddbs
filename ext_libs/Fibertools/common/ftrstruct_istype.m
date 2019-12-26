function [res, errStr, verStr]= ftrstruct_istype(dataIn, fiberType)
%function [res, errStr]= ftrstruct_istype(dataIn, fiberType)
%
%   function determines if dataIn is a valid ftrStruct or fiberStruct. If this is the case,
%   also the type name is returned as the second return value. 
%
%   Type Explanation:
%    ftrStruct is the whole container, containing al fibertracks and all
%        subsets of fibers
%    fiberStruct is a structure containing the information of a fiber
%        subset. In general a ftrStruct contains some fiberstruct in the 'fiber'
%        entry
%
%   dataIn:    arbitary data .. may be a ftrStruct or a fiberStruct
%   fiberType: {'ftrStruct' | 'fiberStruct'} the default value is 'ftrStruct'
%
%   return values
%   res:    1 if arg is an ftrStruct 0 else
%   errStr: If an error occured, errStr is identify the error
%   verStr: String identifying the versrion of the ftrStruct
%
%
% Bjoern W. Kreher
% 11/02
%
% UNIX



res= 0; errStr= ''; typeNameStr= '';, verStr= '';

if ~exist('fiberType') | isempty(fiberType)
    fiberType= 'ftrStruct';
end


if ~isstruct(dataIn)
    errStr= 'Argument is not a struct';
    return
end

nameIndex= fieldnames(dataIn);
if strcmp(fiberType, 'ftrStruct')
    if isfield(dataIn, 'connectCell') & isfield(dataIn, 'curveSegCell') & isfield(dataIn, 'posSegCell') & ...
            isfield(dataIn, 'vox') & isfield(dataIn, 'patient') & isfield(dataIn, 'dtdType') & ...
            isfield(dataIn, 'user') & isfield(dataIn, 'fiber') & isfield(dataIn, 'algoName') & ...
            isfield(dataIn, 'trackParam') & isfield(dataIn, 'trackDate') & isfield(dataIn, 'logData')
        res= 1;
        if isfield(dataIn, 'version') & isfield(dataIn, 'version')
            verStr= dataIn.version;
            return
        else
            verStr= 'V1.0';
            return
        end
    else
        errStr= 'ftrstruct_istype: data is not of the type ftrStruct';
        return
    end
elseif strcmp(fiberType, 'fiberStruct')
    
    nameIdx= find(strcmp(nameIndex, 'name'));
    curveIDIdx=   find(strcmp(nameIndex, 'curveID'));
    roiNameIdx=   find(strcmp(nameIndex, 'roiName')); %: 'cArea_0002'
    userIdx=   find(strcmp(nameIndex, 'user')); %: []
    if ~isempty(nameIdx) & ~isempty(curveIDIdx) & ~isempty(roiNameIdx) & ~isempty(userIdx) & ...
            (isempty(dataIn.name) | isstr(dataIn.name)) & ...
            (isempty(dataIn.curveID) | isnumeric(dataIn.curveID)) & ...
            (isempty(dataIn.roiName) | iscell(dataIn.roiName)) & ...
            ~isempty(userIdx)
        res= 1;
    else
        errStr= 'ftrstruct_istype: data is not of the type fiberStruct';
        return
    end
else
    errStr= sprintf('ftrstruct_istype: ''%s'' is an undefined fiber type', fiberType);
    return
end