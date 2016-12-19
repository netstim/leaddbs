function [res, errStr, oArg1, oArg2]= mrstruct_export(varargin)
%   function [dtdStrcut, errStr, oArg1, oArg2]= mrstruct_export(typeStr, mrType, [, op1[, op2[... [, opN]]]])
%
%   command:    {'SPM02' | '' | 
%
%   'SPM02': input_fName
%
% Bjoern W. Kreher
% 08/05
%
% UNIX

res= []; errStr= ''; oArg1= []; oArg2= [];
maxArg= 10;

if (nargin >= 1) && mrstruct_istype(varargin{1})
    mrStruct= varargin{1};
else
    errStr= strcat(mfilename, ' (error): First argument have to be of the type mrStruct');
    return
end
if (nargin >= 2) && ischar(varargin{2})
    typeStr= varargin{2};
else
    errStr= strcat(mfilename, ' (error): Second argument have to be of the type string');
    return
end

argCell= cell(1, maxArg);
for i= 1:(nargin - 2)
    argCell{i}= varargin{2 + i};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp('SPM02', typeStr)
    headerStrc= argCell{1};
    fName= argCell{2};
    
    [res, errStr]= local_export_SPM02(mrStruct, headerStrc, fName);
else
    errStr= strcat(mfilename, ' (error): Data type ''', typeStr, ''' is not supported yet');
    return    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%% Start:
%%
%%  function [res, errStr]= local_export_SPM02(mrStruct, headerStruc, fName)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_export_SPM02(mrStruct, hdrStrc, fName)
res= []; errStr= '';


sizeAy= mrstruct_query(mrStruct, 'sizeAy');
if length(sizeAy) ~= 3
    errStr= strcat(mfilename, '::local_export_SPM02 (error): Routine supports only volumes yet');
    return
end

hdrStrc.fname= fName;
[res, errStr]= private_writeSPM2(hdrStrc, private_exportSPM2_data(mrStruct.dataAy));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= private_writeSPM2(hdr, data)
res= []; errStr= ''; 

if (exist('spm_vol') ~= 2) || (exist('spm_write_vol') ~= 2)
    errStr= strcat(mfilename, '::private_writeSPM2 (error): Needed SPM02 routines weren''t found');
    return;    
end

try
    res = spm_write_vol(hdr, data);
catch
    errStr= strcat(mfilename, '::private_writeSPM2 (error): SPM02 had trouble by writing a file');
    res= [];
    return;    
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataOut, errStr]= private_exportSPM2_data(dataIn)
dataOut= []; errStr= '';

dataOut= permute(double(dataIn(end:-1:1, :, :)), [2 1 3]);                     % SPM-> mrStruct

