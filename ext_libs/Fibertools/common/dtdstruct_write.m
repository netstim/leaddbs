function [res, errStr]= dtdstruct_write(dtdStruct, fileName)
%
% function [res, errStr]= dtdstruct_write(dtdStruct, fileName)
% 
% writes a dtdStruct to the harddisk
%
% dtdStruct:    should contain a valid dtdStruct
% fileName:     should contain the path and file name. if it is empty or undefined, a file selection dialog will appear 
% 
% return values:
%  res:  if all went fine, res contains the path, where the dtdStruct was saved, otherwise res is empty
%  errStr: if an error occured, errStr contains a string to identify the error
%
%
% Bjoern W. Kreher
% 08/02
%
% UNIX

res= []; errStr= [];

[ok, typeStr, errStr]= dtdstruct_istype(dtdStruct);
if ok && strcmp(typeStr, 'mrStruct')
    dtdStruct= dtdstruct_init('MR', dtdStruct);
end


if ~ok
    errStr= strcat('dtdstruct_write::', errStr);
    return
end

if ~exist('fileName') || isempty(fileName)
    [fileStr,dirStr]=uiputfile('*.mat','select a dtdStruct');
    if isnumeric(fileStr)
        return
    end
    fileName= fullfile(dirStr, fileStr);
end


stNames= fieldnames(dtdStruct);
dtdCell= struct2cell(dtdStruct);

commandStr= sprintf('save(''%s''', fileName);

for i= 1:length(dtdCell)
    commandStr= sprintf('%s, ''%s''', commandStr, stNames{i});
    eval(sprintf('%s= dtdStruct.%s;', stNames{i}, stNames{i}));
end

commandStr= strcat(commandStr, ',''-v7.3'');');
try
    eval(commandStr);
catch
    errStr= strcat('dtdstruct_write (error):', lasterr);
    return
end
res= fileName;


