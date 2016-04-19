function [mrStruct, errStr] = nifti_to_mrstruct(mrType, fNameCell, targetMR, tmpDir)
%
% Import NifTI images into mrstruct
% FORMAT [mrStruct, errStr] = nifti_to_mrstruct(mrType, fNameCell, targetMR, tmpDir)
% ======
% Input arguments
%   mrType    - a valid mrstruct type (all image, volume, series types 
%               are supported, but no spectral types)
%   fNameCell - cell array of filenames. The meaning of cell array 
%               dimensions depends on the mrType:
%                 'image','volume'      - a 1x1 cell array
%                 '(image|volume)Echos' - #Echos-by-1 cell array
%                 'series2D','series3D' - #Timepoints-by-1 cell array
%                 '(series2D|series3D)Echos' - #Echos-by-#Timepoints cell array
%   targetMR - a mrStruct or a nifti file. If defined, fNameCell will
%              resliced to the coordinate system of targrtMR before importing
%   tmpDir   - Directory for temporary files. If not defined the system
%              temporary directory is tried. If this can't be determined,
%              the current directiory is taken
%
% Output arguments
%   mrStruct - mrStruct with all data and position information
%   errStr   - empty, if successful. Otherwise error description.
% All files that have to be imported into a single mrstruct need to have 
% identical voxel sizes, image dimensions and slice orientation. If this 
% is not the case, they need to be resliced within SPM before they can be 
% imported.
% Data are rearranged from NifTI array order to mrStruct array order. A 
% coordinate transformation matrix from mrStruct voxel coordinates to mm 
% coordinates is computed based on the NifTI transformation matrix.
%_______________________________________________________________________
% Bjoern W. Kreher
% 08/05
%
% UNIX
%_______________________________________________________________________
% $Id: nifti_to_mrstruct.m,v 1.10 2011/11/02 11:23:39 glauche Exp $

rev = '$Revision';

mrStruct = []; errStr= '';


% check input parameter
if nargin < 3
    targetMR= [];
end;

if nargin < 4
    try
        tmpDir= tempdir;
        if isempty(tmpDir) || ~exist(tmpDir,'dir')
            tmpDir = pwd;
        end
    catch
        tmpDir = pwd;
    end
end;

if ischar(fNameCell)
    fNameCell= cellstr(fNameCell);
end;

%% reslice files if necessary
if ~isempty(targetMR)
    [fNameCell, errStr, tmpFileCell]= local_reslice(fNameCell, targetMR, tmpDir);
else
    tmpFileCell= {};
end;

%% import files
try
    header = spm_vol(fNameCell);
    spm_check_orientations(cat(1,header{:}));
catch
    errStr = lasterr;
    return;
end;

if isempty(header)
    return
end;
[sizeAy, errStr, templateMRs, mrType]= private_importSPM_header(header, ...
                                                  mrType);
if isempty(templateMRs)
    return
end;

data = reshape(private_importSPM_data(spm_read_vols(cat(1,header{:}))),sizeAy);

mrStruct= mrstruct_init(mrType, data, templateMRs);

%% clean up
for i= 1:length(tmpFileCell)
    delete(tmpFileCell{i}); 
end;

end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, tmpFile]= local_reslice(source, target, tmpDir)
res= []; errStr= ''; tmpFile= {};

% if target is file name of mrStruct --> load mrstruct
if ischar(target)
    [p, f, e]= fileparts(target);
    if strcmp(e, '.mat')
        target= mrstruct_read(target);
        if isempty(target)
            errStr= sprintf('%s::local_reslice(error): could not open mrStruct', mfilename);
            return;
        end;
    end;
end;

% if target is a mrStruct --> convert it to nii
if mrstruct_istype(target)
    fNameTmp= sprintf('%s.nii',tempname(tmpDir));
    [ref, errStr]= mrstruct_to_nifti(target, fNameTmp);
    if isempty(ref)
        return;
    end
    tmpFile{end+1, 1}= fNameTmp;
    target= fNameTmp;
end

% target is now a nifti file (hopefully) reslice bilder

% reslice bilder
resliceFlags.mask= 0;
resliceFlags.meam= 0;
resliceFlags.interp= 1;    % hoffentlich der richtige switch
resliceFlags.which=1;
try 
    spm_reslice({target; source{:}}, resliceFlags);
catch
    errStr= sprintf('spm_reslice: %s', lasterr);
    return;
end

% search temporary files
tmpFile= append_tmpFile(tmpFile, target, 'mean');
tmpFile= append_tmpFile(tmpFile, target, 'r');
tmpFile= append_tmpFile(tmpFile, source, 'r');


% generate resliced filenames
res= cell(length(source), 1);
for i= 1:length(source)
    [p, f, e]= fileparts(source{i});
    res{i}= fullfile(p, strcat('r', f, e));
end;

end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= append_tmpFile(res, fName, preFix)
errStr= '';

if nargin < 3
    preFix= '';
end

if iscell(fName)
    for i= 1:length(fName)
        [temp, errStr]= append_tmpFile({}, fName{i}, preFix);
        res((end+1):(end+length(temp)), 1)= temp(:);
        return
    end
end

[p, f, e]= fileparts(fName);

if exist(fullfile(p, strcat(preFix, f, '.nii')), 'file') == 2
    res{end+1, 1}= fullfile(p, strcat(preFix, f, '.nii'));
end
if exist(fullfile(p, strcat(preFix, f, '.hdr')), 'file') == 2
    res{end+1, 1}= fullfile(p, strcat(preFix, f, '.hdr'));
end
if exist(fullfile(p, strcat(preFix, f, '.img')), 'file') == 2
    res{end+1, 1}= fullfile(p, strcat(preFix, f, '.img'));
end;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataAy, errStr]= private_importSPM_data(data)
dataAy= []; errStr= '';
switch ndims(data)
case 2,
    dataAy= permute(double(data(:, end:-1:1)), [2 1]);                     % SPM-> mrStruct
case 3,
%    dataAy= permute(double(data(:, end:-1:1,end:-1:1)), [2 1 3]);                     % SPM-> mrStruct
    dataAy= permute(double(data(:, end:-1:1,:)), [2 1 3]);                     % SPM-> mrStruct
case 4,
%    dataAy= permute(double(data(:, end:-1:1,end:-1:1,:)), [2 1 3 4]);                     % SPM-> mrStruct
    dataAy= permute(double(data(:, end:-1:1,:,:)), [2 1 3 4]);                     % SPM-> mrStruct
end;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sizeVolAy, errStr, mrStruct, mrType]= private_importSPM_header(header, mrType)
errStr= ''; mrStruct= [];

if exist('spm_imatrix') ~= 2
    errStr= strcat(mfilename, '::private_importSPM_header (error): Needed SPM routines weren''t found');
    mrStruct= [];
    return;    
end

% wegen transponieren bei mrstruct
transMx= [[0 1 0 0]; [1 0 0 0]; [0 0 1 0]; [0 0 0 1]];
trans2Mx= [[1 0 0 0]; [0 -1 0 header{1}.dim(2) + 1]; [0 0 1 0]; [0 0 0 1]];
mrStruct= mrstruct_init;
P= spm_imatrix(header{1}.mat);

mrStruct.vox= abs([P([8 7 9]) 0]);                                          % SPM -> mrStruct

% re-conversion according to spm_dicom_convert
patient_to_tal = diag([-1 -1 1 1]);
%%analyze_to_dicom = [diag([1 -1 1]) [0 (header{1}.dim(2)-1) 0]'; 0 0 0 1]*[eye(4,3) [-1 -1 -1 1]'];
analyze_to_dicom = [diag([1 -1 1]) [0 (header{1}.dim(2) + 1) 0]'; 0 0 0 1];  %%% 28.1.2008 BWK Teil von MVM korregiert

corrMy= diag(ones(1, 4)); corrMy(1:3, 4)= 1;        %%% bei neuer norm muss das rein EDGES_NEW
transMx= transMx*corrMy;                             %%% bei neuer norm muss das rein EDGES_NEW

mrStruct.edges=inv(patient_to_tal)*header{1}.mat*inv(analyze_to_dicom)*transMx;

szhdr = size(header);
szhdr = szhdr(szhdr>1);

if numel(szhdr) > 2
    errStr = strcat(mfilename, ['::private_importSPM_header (error): ' ...
                        'Unsupported number of file array dimensions']);
    mrStruct= [];
    return;
end;

% in SPM, 2D images are just volumes with a z dim of 1
if header{1}.dim(3) ==1
    sizeVolAy = [header{1}.dim([2 1]) szhdr]; % SPM -> mrStruct
    mrType    = strrep(mrType, 'volume', 'image');
    mrType    = strrep(mrType, '3D', '2D');
else
    sizeVolAy= [header{1}.dim([2 1 3]) szhdr];% SPM -> mrStruct
    mrType    = strrep(mrType, 'image', 'volume');
    mrType    = strrep(mrType, '2D', '3D');
end;
end
