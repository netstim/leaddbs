function [res, errStr]= maskstruct_write(mStruct, fileName, diagFlag)
%
% function [res, errStr]= ftrstruct_write(maskStruct, fName)
%
% writes a maskStruct to the harddisk
%
% maskStruct:   should contain a valid maskStruct Ver2
% fileName:     should contain the path and file name. if it is empty or undefined, a file selection dialog will appear 
% 
% return values:
%  res:  if all went fine, res contains the path, where the maskStruct was saved, otherwise res is empty
%  errStr: if an error occured, errStr contains a string to identify the error
%
% Bjoern W. Kreher
% 12/07
%
% UNIX

res= []; errStr= [];

%% check maskStruct
if ~maskstruct_istype(mStruct)
    errStr= 'maskstruct_write: First argumenst is not from type maskStruct Ver2';
    return;
end

%% check resp. determine file name
if ~exist('fileName') || isempty(fileName)
    [fileStr,dirStr]=uiputfile('*.mat','select a maskStruct');
    if isnumeric(fileStr)
        return
    end
    fileName= fullfile(dirStr, fileStr);
elseif exist('diagFlag') && strcmp(diagFlag, 'yes')
    [fileStr,dirStr]=uiputfile(fileName,'select a maskStruct');
    if isnumeric(fileStr)
        return
    end
    fileName= fullfile(dirStr, fileStr);    
end

if ~ischar(fileName)
    errStr= 'maskstruct_write: Filename from type string';
    return;
end

%% save data
maskCell= mStruct.maskCell;
maskNamesCell= mStruct.maskNamesCell;
sizeAy= mStruct.sizeAy;
mrsProp= mStruct.mrsProp;
user= mStruct.user;
version= mStruct.version;

save(fileName, 'maskCell', 'maskNamesCell', 'sizeAy', 'mrsProp', 'user', 'version');

res= fileName;