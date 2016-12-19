function [res, errStr, oArg1, oArg2]= mrstruct_import(varargin)
%   function [dtdStrcut, errStr, oArg1, oArg2]= mrstruct_import(typeStr, mrType, [, op1[, op2[... [, opN]]]])
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


if (nargin >= 1) && ischar(varargin{1})
    typeStr= varargin{1};
else
    errStr= strcat(mfilename, ' (error): Second argument have to be of the type string');
    return
end

argCell= cell(1, maxArg);
for i= 1:(nargin - 1)
    argCell{i}= varargin{1 + i};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp('SPM02', typeStr)
    [res, errStr]= local_import_SPM02(argCell{1}, argCell{2});
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
%%  function [mrStruct, errStr]= local_import_n_SPM2(header, data)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_import_SPM02(mrType, fNameCell)
res= []; errStr= '';

if exist('spm') ~= 2
    errStr= strcat(mfilename, '::local_import_SPM02 (error): Could not find SPM');
    return;    
end
if ~strcmp('SPM2', spm('ver'))
    errStr= strcat(mfilename, '::local_import_SPM02 (error): Found wrong SPM version');
    return;    
end

if ischar(fNameCell)
    fNameCell= {fNameCell};
end
[header, errStr, data]= private_openSPM2(fNameCell{1});
if isempty(header)
    return
end
[sizeVolAy, errStr, templateMRs]= private_importSPM2_header(header);
if isempty(templateMRs)
    return
end

if numel(fNameCell) == 1 % if volume
    res= mrstruct_init(mrType, private_importSPM2_data(data), templateMRs);
    return
end

% bestimme sizeAy
sizeAy= [sizeVolAy size(fNameCell)];
volVoxNo= prod(sizeVolAy);

mrStruct= mrstruct_init(mrType, zeros(sizeAy), templateMRs);
mrStruct.dataAy(1:volVoxNo)=  private_importSPM2_data(data);

for i= 2:numel(fNameCell)
    [header, errStr, data]= private_openSPM2(fNameCell{i});
    if isempty(header)
        return
    end
    sizeAyNew= private_importSPM2_header(header);
    if ~isequal(sizeVolAy, sizeAyNew)
        errStr= strcat(mfilename, '::local_import_SPM02 (error): SPM file ''', fNameCell{i}, ''' has different size');
        return
    end
    mrStruct.dataAy(((i - 1)*volVoxNo + 1):(i*volVoxNo))= private_importSPM2_data(data);
end
 
res= mrStruct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hdr, errStr, data]= private_openSPM2(fNameIn)
hdr= []; errStr= ''; data= [];

[pStr, fStr, extStr]= fileparts(fNameIn);
fName= fullfile(pStr, strcat(fStr, '.img'));

if ~exist(fName, 'file')
    errStr= strcat(mfilename, '::private_openSPM2 (error): File ''', fName, ''' not found');
    return;    
end

if (exist('spm_vol') ~= 2) || (exist('spm_read_vols') ~= 2)
    errStr= strcat(mfilename, '::private_openSPM2 (error): Needed SPM02 routines weren''t found');
    return;    
end

try
    hdr = spm_vol(fName);
catch
    errStr= strcat(mfilename, '::private_openSPM2 (error): SPM02 had trouble with file ''', fName, '''');
    hdr= [];
    return;    
end
    
try
    data= spm_read_vols(hdr);
catch
    errStr= strcat(mfilename, '::private_openSPM2 (error): SPM02 had trouble with file ''', fName, '''(header seemed to be ok)');
    hdr= [];    data= [];
    return;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataAy, errStr]= private_importSPM2_data(data)
dataAy= []; errStr= '';

dataAy= permute(double(data(:, end:-1:1, :)), [2 1 3]);                     % SPM-> mrStruct
%dataAy= double(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sizeVolAy, errStr, mrStruct]= private_importSPM2_header(header)
errStr= ''; mrStruct= [];

if exist('spm_imatrix') ~= 2
    errStr= strcat(mfilename, '::private_importSPM2_header (error): Needed SPM02 routines weren''t found');
    mrStruct= [];
    return;    
end

if nargout == 3
    % wegen transponieen bei mrstruct
    transMx= [[0 1 0 0]; [1 0 0 0]; [0 0 1 0]; [0 0 0 1]];
    trans2Mx= [[1 0 0 0]; [0 -1 0 header.dim(2) + 1]; [0 0 1 0]; [0 0 0 1]];
    mrStruct= mrstruct_init;
    P= spm_imatrix(header.mat);
    mrStruct.user.spm02_Header= header;
    
    mrStruct.user.hMatrix= header.mat*trans2Mx*transMx;                         % SPM -> mrStruct
    mrStruct.edges= mrStruct.user.hMatrix;
    mrStruct.vox= abs([P([8 7 9]) 0]);                                          % SPM -> mrStruct
%    mrStruct.user.hMatrix= header.mat;
%    mrStruct.vox= abs([P([7 8 9]) 0]);
end

sizeVolAy= header.dim([2 1 3]);                                             % SPM -> mrStruct
%sizeVolAy= header.dim([1 2 3]);
