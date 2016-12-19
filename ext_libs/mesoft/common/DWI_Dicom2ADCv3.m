%   DWI_Dicom2ADCv2 
%
%   converts DW Dicom Images in mr-sructure and calculates DTI (DTD structure)
%   ika 031030
% 
%   y = DWI_Dicom2ADCv3(Scans, N_slices,  b_factor, tresholdAbs, pathname1, filename1)
%   input parameters
%               Scans               scan(measurement) number , can be a vector, for example [17, 19, 21:22] 
%                                   it is assumesd that all scans are identical., i.e. averages over the scans list
%               N_slices            number of slices in each serie
%               b_factor            b-factor, can be also vestor of several differnt b-factors, as for example [0, 200, 555, 999, 1200] 
%
%   output parameters
%               ExecutionTime       time of all calculation
%               the DTD structure will be written in the "patient home directory"
%
%
%   example  ExecutionTime = DWI_Dicom2ADCv2( (17:22), 52, [0,750,1000]))
%
%   Wriiten by:
%   ------------------------------------------------------------
%   Dr.Kamil Il'yasov
%   Section of Medical Physics, Dept. of Diagnostic Radiology
%   University Medical Center of Freiburg
%   Hugstetter Str. 55, D-79106 Freiburg, Germany
%   Tel.: +49 761 270 7392
%   Fax : +49 761 270 3831
%   E_mail: Kamil.Ilyasov@uniklinik-freiburg.de
%   ------------------------------------------------------------
%
% ika031116     added treshold, all  pixels on the 1st B0 image having intensity below this treshold will be excluded from calculation and set to 0
% ika 050126    added option for scripting and export in ROI_tool compartible format 



function y = DWI_Dicom2ADCv3(Scans, N_slices,  b_factor, tresholdAbs, pathname1, filename1)


tic

if (nargin <3) 
    'Not enouth input parameters, program can not be executed '
end;

if (nargin <4 ) 
tresholdAbs =40;
end

if (nargin <6 ) 
[filename1, pathname1] = uigetfile('*.dcm', 'Select a file corresponding to the first DWI scan');
end


N_DWIperSlice =length(b_factor);

index1 = find(pathname1 =='\' | pathname1 =='\');
FilesRootDirectory =pathname1(1:(index1(end-1))); % % 'E:\zDataTemp\albino7'

index2 = find(filename1 =='s' );
PatientName = filename1(1:(index2(end)-2) ); % 'Cerny_Esther_20030926' %'Nedeljko_Antinic_20030918'  edynak_Andreas_20031022_s017_i0002.dcm

%

for i=1:length(Scans)
    ScanString = num2str(Scans(i));
    if length(ScanString) ==1, ScanString=['000',ScanString]; end
    if length(ScanString) ==2, ScanString=['00',ScanString]; end
    if length(ScanString) ==3, ScanString=['0',ScanString]; end
    FileName= [ PatientName,'_s',ScanString,'_i0001.dcm'];
    
    CurrentScanDirectory =[FilesRootDirectory,  ScanString, '\'];
    % res = DicomReadSerie(N_slices, 13,FileName , CurrentScanDirectory,FilesRootDirectory );
    res = DicomReadSerieV2(N_slices, N_DWIperSlice,FileName , CurrentScanDirectory,FilesRootDirectory );
end
% 

tmp1=[];
for i=1:length(Scans)
    ScanString = num2str(Scans(i));
    if length(ScanString) ==1, ScanString=['000',ScanString]; end
    if length(ScanString) ==2, ScanString=['00',ScanString]; end
    if length(ScanString) ==3, ScanString=['0',ScanString]; end
    text2=[' load ', FilesRootDirectory, 'DTI_',  PatientName, '_s',ScanString];
    eval(text2);
    if i==1
        tmp1 = DTI.dataAy;
    else
        tmp1 = tmp1+ DTI.dataAy;
    end
end

DTI.dataAy = (1/length(Scans))*tmp1;
text3 =['save ', FilesRootDirectory,'DTI_',PatientName, '_s',num2str(Scans(1)),'_', num2str(Scans(length(Scans)))];
eval(text3);

NameOfOutputStructure = [PatientName, '_s',num2str(Scans(1)),'_', num2str(Scans(length(Scans)))]
DE_iso(NameOfOutputStructure,FilesRootDirectory, DTI, N_slices, tresholdAbs, b_factor)

y =toc; % calculation time
