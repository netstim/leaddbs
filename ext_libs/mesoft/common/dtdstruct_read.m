function [res, errStr, fName]= dtdstruct_read(fName)
%
%function [res, errStr, fName]= dtdstruct_read(fName)
%
%   Opens a dtdStruct from the filesystem
%
% fName:     should contain the path and file name. if it is empty or undefined, a file selection dialog will appear 
% 
% return values:
%  res:  contains the opend dtdStruct. If an error occurs res is empty
%  errStr: if an error occured, errStr contains a string to identify the error
%  fName:  if all went fine, fName contains the path of the opend dtdstruct
%
%
% Bjoern W. Kreher
% 08/02
%
% UNIX

errStr= '';

if nargin == 0 
    [fileStr,dirStr]=uigetfile('*.mat','Load a dtdStruct');
    if fileStr == 0
        fName= [];
    else
        fName= fullfile(dirStr, fileStr);
    end
end
if exist(fName) == 2
    res= open(fName);
    [ok, typeStr, errStr]= dtdstruct_istype(res);
    if ok && strcmp(typeStr, 'mrStruct')
        res= dtdstruct_init('MR', res);
    end

    if strcmp(typeStr, 'soft-link dtdStruct')
        if strcmp(res.PRIV_FORMAT, 'SPM02')
            [res, errStr]= dtdstruct_import('priv_SPM02_SL', fName);
            return;
        else
            errStr= strcat(mfilename,' (error): The format ''', res.PRIV_FORMAT, ''' is not supported as softlink yet');  res= [];
            return;
        end
    end
    if ~dtdstruct_istype(res)
        res= mrstruct_read(fName);
        if mrstruct_istype(res) 
            [res, errStr]= dtdstruct_init('MR', res);
        else
            res= [];
            errStr= 'dtdstruct_read: file is not from type dtdstrtuct';
        end
    end
else
    res= [];
    errStr= 'file not found';
end
