function [res, errStr]= fiberstruct_write(fiberStruct, fileName)
%
% function [res, errStr]= fiberstruct_write(fiberStruct, fileName)
%
%
%  Bjoern W. Kreher
%  08/02
%
%  UNIX
%

res= []; errStr= [];

if ~exist('fileName')
    [fileStr,dirStr]=uiputfile('*.mat','select a fiberStruct');
    if isnumeric(fileStr)
        return
    end
    fileName= strcat(dirStr, '/', fileStr);
end

curvesData= fiberStruct.curvesData;
save(fileName, 'curvesData');

res= fileName;