function [res, errStr, fName]= ftrstruct_read(fNameIn)
%
%function [res, errStr, fName]= ftrstruct_read(fName)
%
%   Opens a ftrStruct from the filesystem
%
% fName:     should contain the path and file name. if it is empty or undefined, a file selection dialog will appear 
% 
% return values:
%  res:  contains the opend ftrStruct. If an error occurs res is empty
%  errStr: if an error occured, errStr contains a string to identify the error
%  fName:  if all went fine, fName contains the path of the opend ftrStruct
%
%
% Bjoern W. Kreher
% 11/02
%
% UNIX

res= []; errStr= ''; fName= '';

if (nargin == 0) || isempty(fNameIn) || (exist(fNameIn, 'dir') == 7)
    if (nargin > 0) & (exist(fNameIn, 'dir') == 7)
        curDir= pwd;
        cd(fNameIn);
    end
    [fileStr,dirStr]=uigetfile('*.mat','Load a ftrStruct');
    if (nargin > 0) & (exist(fNameIn, 'dir') == 7)
        cd(curDir);
    end

    if fileStr == 0
        fNameIn= [];
    else
        fNameIn= fullfile(dirStr, fileStr);
    end
end
fName= fNameIn;
if exist(fNameIn, 'file') == 2
    res= open(fNameIn);
    % work around for the old ftrStruct format
    if isstruct(res)
        fNames= fieldnames(res);
        if isempty(find(strcmp(fNames, 'algoName')))
            res.algoName= '';
            res.trackParam= [];
            res.trackDate= '';
            res.logData= [];
        end
    end
else
    res= [];
    errStr= 'ftrStruct_read: file not found';
    return
end

if ~ftrstruct_istype(res)
    names= fieldnames(res);
    if ~isempty(find(strcmp(names, 'curvesData'))) & ~isempty(find(strcmp(names, 'voxelData')))
        ftrStruct= ftrstruct_init;
        ftrStruct.dtdType= 'DTD';
        ftrStruct.posSegCell= res.voxelData;
        ftrStruct.curveSegCell= res.curvesData;
        for i= 1:length(ftrStruct.curveSegCell)
            ftrStruct.connectCell{i}= [i];
        end
        res= ftrStruct;
    else
        res= [];
        errStr= 'ftrStruct_read: File was not from type ftrStruct';
    end
end

