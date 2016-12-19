% DWI_Dicom2DTD 
% version 3.01
% Converts DW Dicom Images in mr-sructure and calculates DTI (DTD structure)
% ika 031030
%y = DWI_Dicom2DTDv3(Scans, N_slices, DT_Encoding_Schema, b_factor, tresholdAbs, pathname1, filename1) 
% input parameters
%               Scans               scan(measurement) number , can be a vector, for example [17, 19, 21:22] 
%               N_slices            number of slices in each serie
%               DT_Encoding_Schema  'DE12' or 'DE06', i.e. currently supports only standard Siemens diffusion encoding direction 
%       `       b_factor            b-factor, if only 3 first parameters are given will be set to 1000
%               tresholdAbs         all pixels in the 1st b0-image below
%                                   tresholdAbs will be ignored, if only 3 or 4 first parameters are given will be set to 40
%                       if less than 7 input parameters are given, the following parameters will be set per dialog
%               pathname1           directory with DWI dicom images ="patient home directory"
%               filename1           name of one of dimom images in the
%                                   serie (numiration of images should start from 1, all images should be available in this directory)
% output parameters
%               ExecutionTime       time of all calculation
%               the DTD structure will be written in the "patient home directory"
%
%
% example1  ExecutionTime = DWI_Dicom2DTDv4( (17:22), 52, 'DE12') ; other
%                                   parameters set to   b_factor = 1000; tresholdAbs =40; pathname1 &
%                                   filename1 will be selected per dialog
% example2 - all parameters set explicitly 
%           ExecutionTimeScans = DWI_Dicom2DTDv3( [28:33] , 52, 'DE12',  1000, 40, 'D:\MyDirectory', 'HALI_EMKA_20040607_s0083_i0001.dcm ') , 
% example3 - uses exact b_tensor 
%           ExecutionTimeScans = DWI_Dicom2DTDv3( [28:33] , 52, 'DE12',  b_tensor) , 
%                                   in this case b_tensor is a matrix of size 3:3:13, i.e. 3x3 b-matrixis for b0
%                                   and 12DE directions
%
%Writen by:
%   ------------------------------------------------------------
%   Dr.Kamil Il'yasov
%   Section of Medical Physics, Dept. of Diagnostic Radiology
%   University Medical Center of Freiburg
%   Hugstetter Str. 55, D-79106 Freiburg, Germany
%   Tel.: +49 761 270 7392
%   Fax : +49 761 270 3831
%   E_mail: Kamil.Ilyasov@uniklinik-freiburg.de
%   ------------------------------------------------------------
% Changes by ika: 
% 031116 added threshold, all pixels on the 1st B0 image having intensity below this threshold will be excluded from calculation and set to 0
% version 3.0 includes calculation for exact b-factor
% ika 050620 rotated grads dir with slice
% 
% Corrected by Susanne Schnell, March 2007



function y = DWI_Dicom2DTDv4(Scans, N_slices, DT_Encoding_Schema, b_factor, tresholdAbs, pathname1, filename1, vHd)

y= [];
tic
if (nargin <3 ) 
    DT_Encoding_Schema = 'DE12';
    b_factor = 1000;
end;
if (nargin <4 ) 
    b_factor = 1000;
end;

if (nargin <2) 
    'Not enouth input parameters, program can not be executed '
end;

if (nargin <5 ) 
tresholdAbs =40;
end

if (strcmp(DT_Encoding_Schema, 'DE12'))
    N_DWIperSlice = 13;
elseif (strcmp(DT_Encoding_Schema,'DE06'))
    N_DWIperSlice = 7;
else
    disp('unknown diffusion ensoding schema, only Siemens 6 & 12 diffusion encoding directions supported')
    return
end

if ~exist('pathname1') || ~exist('filename1') || isempty(pathname1) || isempty(filename1)
    [filename1, pathname1] = uigetfile('*.dcm', 'Select a file corresponding to the first DWI scan');
    drawnow;
end

if ~exist('vHd') || ~ishandle(vHd)
    vHd= [];
end

% a bit more complex for arbitrary b-matrix od multiple b-factors
b_factorSize = size(b_factor);
if (length(b_factor)>1 )
    if (length(b_factorSize)==3 && b_factorSize(1)==3 && b_factorSize(2)==3) 
        N_DWIperSlice = length(b_factor);
    elseif(length(b_factorSize)==2 && b_factorSize(1)==1 ) 
        if b_factor(1)==0
            N_DWIperSlice = 1+ (N_DWIperSlice-1)*(length(b_factor) -1)
        else
            N_DWIperSlice = 1+ (N_DWIperSlice-1)*(length(b_factor) -1)
        end
    else
        'b-factor has unsupported format'
        ret = 0;
        return
    end
end;


index1 = find(pathname1 =='\' | pathname1 =='/');
FilesRootDirectory =pathname1(1:(index1(end-1)));

index2 = find(filename1 =='s' );
PatientName = filename1(1:(index2(end)-2) );

index3 = find(filename1 =='i' );
index4 = find(filename1 =='.' );
% check how long is the field for scan number
ScanStringLength = index3(end) - index2(end) - 2;
ImageStringLength = index4(end) - index3(end) - 1;


for i=1:length(Scans)
    ScanString= sprintf(sprintf('%%0%dd', ScanStringLength), Scans(i));
    ImageString= sprintf(sprintf('%%0%dd', ImageStringLength), 1);
    FileName= [ PatientName,'_s',ScanString,'_i', ImageString, '.dcm']; 
    CurrentScanDirectory =[FilesRootDirectory,  ScanString, filesep];
    if ~exist(fullfile(CurrentScanDirectory, FileName), 'file') 
        CurrentScanDirectory =pathname1;
        if ~exist(fullfile(CurrentScanDirectory, FileName), 'file')
            errStr= sprintf('%s (error): File %s does not exist', mfilename, fullfile(CurrentScanDirectory, FileName));
            if isempty(vHd)
                disp(errStr);
            else
                set(vHd, 'String', errStr);
                drawnow;
            end
            return
        end
    end
    msgStr= sprintf('Prepare Data %d/%d', i, length(Scans));
    if isempty(vHd)
        disp(msgStr);
    else
        set(vHd, 'String', msgStr);
        drawnow;
    end   
    [DE_dirs,newFile(i,:)] = DicomReadSerieV2(N_slices, N_DWIperSlice,FileName , CurrentScanDirectory,FilesRootDirectory );
    DE_dirs = DE_dirs';
end

tmp1=[];
for i=1:length(Scans)
    DTI = mrstruct_read(newFile(i,:));
    if i==1       
        tmp1 = DTI.dataAy;
    else
        tmp1 = tmp1+ DTI.dataAy;
    end
end

if Scans > 1
    DTI.dataAy = (1/length(Scans))*tmp1;
    text3 =['save(''', FilesRootDirectory,'DTI_',PatientName,'_s',num2str(Scans(1)),'_', num2str(Scans(length(Scans))), ''');'];
    eval(text3);
    NameOfOutputStructure = [PatientName, '_s',num2str(Scans(1)),'_', num2str(Scans(length(Scans)))];
else
    NameOfOutputStructure = [PatientName, '_s',num2str(Scans(1))];
end


msgStr= sprintf('Starting tensor calculation ...');
if isempty(vHd)
    disp(msgStr);
else
    set(vHd, 'String', msgStr);
    drawnow;
end

%DE_schema12_b0_1000(NameOfOutputStructure,FilesRootDirectory, DTI, N_slices, tresholdAbs, DT_Encoding_Schema,b_factor)
DE_schema12_b0_1000(NameOfOutputStructure,FilesRootDirectory, DTI, N_slices, tresholdAbs, DE_dirs,b_factor)

msgStr= sprintf('Done');
if isempty(vHd)
    disp(msgStr);
else
    set(vHd, 'String', msgStr);
    drawnow;
end
y =toc; % calculation time
