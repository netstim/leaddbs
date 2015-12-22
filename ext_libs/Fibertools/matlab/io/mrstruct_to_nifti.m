function [res, errStr]= mrstruct_to_nifti(mrStruct, fName, dataType)
% Convert mrStruct image/volume/series to NifTI compatible files
% FORMAT [res, errStr]= mrstruct_to_nifti(mrStruct, fName, dataType)
% ======
% This routine tries to convert mrStruct data to NifTI compatible files 
% using SPM functions spm_type, spm_platform, spm_write_vol. Image data 
% will be saved in user given format (dataType). Default format is int16. 
% If available, a NifTI compatible coordinate transformation matrix will be
% computed from mrStruct.edges. If this is not possible, axial orientation
% is assumed and voxel sizes are read from mrStruct.vox. If no voxel sizes
% are given, a default of 1x1x1mm will be assumed. Image data are
% rearranged in SPM/Analyze format and slices are re-ordered depending on
% the handedness of the voxel-to-world coordinate mapping.
%
% Input arguments
% mrStruct - mrstruct to be converted
% fName    - location and filename for created file(s). If mrStruct 
%            contains series or multi-echo data, separate volumes will be 
%            created for each echo or series image/volume by appending a 
%            running series/echo number to the file name.
% dataType - data type for output ('uint8','int16','int32','float32',
%            'float64','int8','uint16','uint32')
% Output arguments
% res      - a cell array containing volume handles (see spm_vol). The 
%            array is shaped according to the series/echo dimensions of 
%            the mrStruct data array.
% errStr   - empty, if successful operation. Otherwise, an error message.
%_______________________________________________________________________
% Bjoern W. Kreher
% 08/05
%
% UNIX
%_______________________________________________________________________
% $Id: mrstruct_to_nifti.m,v 1.14 2013/07/26 11:16:11 reisertm Exp $

rev = '$Revision';

res= []; errStr= '';

if exist('spm', 'file') 
    spmVer= spm('ver');
else
    spmVer= '';
end;

if ~any(strcmp(spmVer, {'SPM5', 'SPM8', 'SPM8b', 'SPM12a', 'SPM12b', 'SPM12'}))
    errStr= 'nifti_to_mrstruct: no compatible SPM version was found';
    return;
end;

if ~mrstruct_istype(mrStruct)
    errStr= strcat(mfilename, ' (error): First argument must be of type mrStruct');
    return
end

dataTypeOut = 'int16';
if nargin == 3
    if isnan(spm_type(dataType))
       warning(['mrstruct_to_nifti: unknown data type ' dataType ' ... exporting uint16']);
    else
        dataTypeOut = dataType;
    end
end

% get data dimensionality
dim = mrstruct_query(mrStruct, 'sizeAy');
imtype = mrstruct_query(mrStruct, 'dataType');
switch imtype
case 'image'
    hdrStrc.dim = [dim([2 1]) 1];
    nser = 1;
    nech = 1;
case 'imageEchos'
    hdrStrc.dim = [dim([2 1]) 1];
    nser = 1;
    nech = dim(3);
case 'series2D'
    hdrStrc.dim = [dim([2 1]) 1];
    nser = dim(3);
    nech = 1;
case 'series2DEchos'
    hdrStrc.dim = [dim([2 1]) 1];
    nser = dim(4);
    nech = dim(3);
case 'volume'
    hdrStrc.dim = dim([2 1 3]);
    nser = 1;
    nech = 1;
case 'volumeEchos'
    hdrStrc.dim = dim([2 1 3]);
    nser = 1;
    nech = dim(4);
case 'series3D'
    hdrStrc.dim = dim([2 1 3]);
    nser = dim(4);
    nech = 1;
case 'series3DEchos'
    hdrStrc.dim = dim([2 1 3]);
    nser = dim(5);
    nech = dim(4);
otherwise
    errStr= strcat(mfilename,...
                   ' (error): mrStruct datatype not supported: ', imtype);
    return
end;

edges = mrstruct_query(mrStruct,'edges');
if isempty(edges),
    warning(['%s: No information about spatial orientation '...
                 'found in mrStruct, assuming default orientation.'], mfilename);
    pos=[eye(3) -hdrStrc.dim(:)/2; 0 0 0 1];
    vox = mrstruct_query(mrStruct,'vox');
    if isempty(vox),
        hdrStrc.mat=pos;
        warning(['%s: No information about voxel size '...
                 'found in mrStruct, assuming default voxel size.'], mfilename);
    else
        hdrStrc.mat=diag([vox([2 1]) sum(vox(3:end)) 1])*pos;
    end;
    flipz = 0;
else
    % wegen transponieren bei mrstruct
    transMx= [[0 1 0 0]; [1 0 0 0]; [0 0 1 0]; [0 0 0 1]];
    % re-conversion according to spm_dicom_convert
    patient_to_tal = diag([-1 -1 1 1]);
%%%    analyze_to_dicom = [diag([1 -1 1]) [0 (hdrStrc.dim(2)-1) 0]'; 0 0 0 1]* ...
%%%        [eye(4,3) [-1 -1 -1 1]'];
    analyze_to_dicom = [diag([1 -1 1]) [0 (hdrStrc.dim(2) + 1) 0]'; 0 0 0 1];  %%% 28.1.2008 BWK Teil von MVM korregiert
    
    corrMy= diag(ones(1, 4)); corrMy(1:3, 4)= -1;        %%% bei neuer norm muss das rein EDGES_NEW

%     % SPM will reorder slices during conversion into ascending order,
%     % mrStruct does not have a convention - see if we need to reorder
%     % slices by checking handedness of coordinate system
%     if det(mrStruct.edges(1:3,1:3)) > 0,
%         flipz = 1;
% %        transdim3 = [diag([1 1 -1]) [0 0 (hdrStrc.dim(3)-1)]'; 0 0 0 1]; 
%         transdim3 = [diag([1 1 -1]) [0 0 (hdrStrc.dim(3) + 1)]'; 0 0 0 1]; %%% von BWK 24.1.2008
%     else
%         flipz = 0;
%         transdim3 = eye(4);
%     end;
    
    

% disabled above code due to several problems, however,cd has to 
% handled with caution, not sure whether above code is really necessary
% (M.Reisert Jul2013)

    flipz = 0;
    transdim3 = eye(4);
        
    
    
%    hdrStrc.mat=patient_to_tal*mrStruct.edges*inv(transMx)*inv(transdim3)*analyze_to_dicom;
    hdrStrc.mat=patient_to_tal*mrStruct.edges*corrMy*inv(transMx)*inv(transdim3)*analyze_to_dicom;      %%% bei neuer norm muss das rein EDGES_NEW
end;


[p fileformat ext v] = spm_fileparts(fName); % get filename and ext part of fname

if ~(strcmp(ext,'.nii')||strcmp(ext,'.img')),
    ext = '.img';
end;

if nech > 1,
    fileformat = sprintf('%s-e%%0.%dd', fileformat, floor(log10(nech))+1);
end;
if nser > 1,
    fileformat = sprintf('%s-t%%0.%dd', fileformat, floor(log10(nser))+1);
end;

hdrStrc.dt = [spm_type(dataTypeOut), spm_platform('bigend')];
hdrStrc.pinfo = [Inf Inf 0]';
res = cell(nech, nser);
for ec = 1:nech
    for se = 1:nser
        switch imtype
        case {'image','imageEchos','series2DEchos'}
            hdrStrc.fname= fullfile(p, [sprintf(fileformat, ec, se) ...
                                ext]);
            res{ec,se} = spm_write_vol(hdrStrc, ...
                                       private_exportSPM_data(mrStruct.dataAy(:,:,ec,se),...
                                                              flipz));
        case 'series2D',
            hdrStrc.fname= fullfile(p, [sprintf(fileformat, se) ...
                                ext]);
            res{se} = spm_write_vol(hdrStrc, ...
                                    private_exportSPM_data(mrStruct.dataAy(:,:,se),...
                                                           flipz));
        case {'volume','volumeEchos','series3DEchos'}
            hdrStrc.fname= fullfile(p, [sprintf(fileformat, ec, se) ...
                                ext]);
            res{ec,se} = spm_write_vol(hdrStrc, ...
                                       private_exportSPM_data(mrStruct.dataAy(:,:,:,ec,se),...
                                                              flipz));
        case 'series3D',
            hdrStrc.fname= fullfile(p, [sprintf(fileformat, se) ...
                                ext]);
            res{se} = spm_write_vol(hdrStrc, ...
                                    private_exportSPM_data(mrStruct.dataAy(:,:,:,se),...
                                                              flipz));
        end;
    end;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataOut, errStr]= private_exportSPM_data(dataIn,flipz)
dataOut = []; errStr= '';
if flipz,
    dataOut= permute(double(dataIn(end:-1:1, :, end:-1:1)), [2 1 3]); % SPM-> mrStruct
else
    dataOut= permute(double(dataIn(end:-1:1, :, :)), [2 1 3]);
end;
end