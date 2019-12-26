function [res, errStr]= ftrstruct_write(ftrStruct, fileName)
%
% function [res, errStr]= ftrstruct_write(ftrStruct, fileName)
%
% writes a ftrStruct to the harddisk
%
% ftrStruct:    should contain a valid ftrStruct
% fileName:     should contain the path and file name. if it is empty or undefined, a file selection dialog will appear 
% 
% return values:
%  res:  if all went fine, res contains the path, where the ftrStruct was saved, otherwise res is empty
%  errStr: if an error occured, errStr contains a string to identify the error
%
% Bjoern W. Kreher
% 11/02
%
% UNIX

res= []; errStr= [];

if ~ftrstruct_istype(ftrStruct)
    errStr= 'ftrstruct_write: First argumenst is not from type ftrstruct';
    return;
end

if ~exist('fileName') | isempty(fileName)
    [fileStr,dirStr]=uiputfile('*.mat','select a ftrStruct');
    if isnumeric(fileStr)
        return
    end
    fileName= fullfile(dirStr, fileStr);
end

if ~isstr(fileName)
    errStr= 'ftrstruct_write: Filename from type string';
    return;
end

connectCell= ftrStruct.connectCell;
curveSegCell= ftrStruct.curveSegCell;
posSegCell= ftrStruct.posSegCell;

vox= ftrStruct.vox;
patient= ftrStruct.patient;
user= ftrStruct.user;
dtdType= ftrStruct.dtdType;

fiber= ftrStruct.fiber;
algoName= ftrStruct.algoName;
trackParam= ftrStruct.trackParam;
trackDate= ftrStruct.trackDate;
logData= ftrStruct.logData;

if strcmp(ftrstruct_query(ftrStruct, 'getVer'), 'V1.0')
    save(fileName, 'connectCell', 'curveSegCell', 'posSegCell', 'fiber', 'dtdType', 'vox', ...
        'patient', 'user', 'algoName', 'trackParam', 'trackDate', 'logData');
elseif strmatch('V1.1', ftrstruct_query(ftrStruct, 'getVer')) == 1
    hMatrix= ftrStruct.hMatrix;
    version= ftrStruct.version;
    save(fileName, 'connectCell', 'curveSegCell', 'posSegCell', 'fiber', 'dtdType', 'vox', ...
        'patient', 'user', 'algoName', 'trackParam', 'trackDate', 'logData', 'hMatrix', 'version');    
end
%save(fileName, 'user', 'patient', 'vox', 'posSegCell', 'curveSegCell', 'connectCell');

res= fileName;