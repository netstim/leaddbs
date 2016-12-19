%DicomReadSerieV2 - reads serie of dicom images and saves in mrstruct
% function res = DicomReadSerie(NSlices, N_de_steps, filename, pathname)
% example res = DicomReadSerie(16, 60)
% ika 031030 - added possibility to write result in arbitrary directory
% ika 030227
% ika030319: bug correction & added check for number of arguments
%
%   Writen by:
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
%   Corrected by Susanne Schnell, March 2007


function [DE_dirs,newFile] = DicomReadSerieV2(NSlices, N_de_steps, filename, pathname, outputPathname)

%check for number of arguments
if (nargin==2) 
    [filename, pathname] = uigetfile('*.dcm', 'Select a file corresponding to the proper scan # ');
end;
if (nargin==3) 
    pathname = pwd; pathname =[pathname, '\'];
end
if (nargin<2) 
    error('not enough input variables, specify at least N-slices and N-repetitions'); 
end;

if (nargin<5) 
    outputPathname = pathname;
end;
% end of errors and tests


[DTI, infoMx]=dicom_read_singlefile([pathname ,filename], 1);

[p, fName, ext]= fileparts(filename);

idx= strfind(fName, '_i') + 2;
idx= idx(end);
digitNo= length(fName) - idx + 1;
for i = 1 : length(infoMx.Private_0029_1010) 
    if strcmp(infoMx.Private_0029_1010(i).name,'NumberOfImagesInMosaic')
        strucNumberOfImagesInMosaic=infoMx.Private_0029_1010(i);
    end
end
if isempty(strucNumberOfImagesInMosaic.item) 
    NumberOfImagesInMosaic = 1; 
    NumberOfImagesInRaw =1; % by mosaic is always NxN pictures
    sizeIm =size(DTI.dataAy);
    mosaic = 0;
else
    NumberOfImagesInMosaic = str2num(strucNumberOfImagesInMosaic.item(1).val);
    tmp = fix(sqrt(NumberOfImagesInMosaic));
    while tmp^2 < NumberOfImagesInMosaic
        tmp = tmp +1;
    end
    if tmp^2 >= NumberOfImagesInMosaic
        NumberOfImagesInRaw =tmp;
    end
    if NumberOfImagesInMosaic < NSlices || NumberOfImagesInMosaic > NSlices
        NSlices = NumberOfImagesInMosaic; % slices given in dti_calulation gui does not agree with existing in mosaic images, correction!
    end
    sizeImOrig =size(DTI.dataAy);
    sizeImOrig(1) = sizeImOrig(1)/NumberOfImagesInRaw ;sizeImOrig(2)= sizeImOrig(2)/NumberOfImagesInRaw ; % size of each subscell in mosaic
    sizeIm =sizeImOrig;
    mosaic = 1;
end

DTI.dataAy =[];
DTI.dataAy =zeros( sizeIm(1),sizeIm(2) ,NSlices, N_de_steps);
if mosaic == 0
    for Slice= 1:NSlices
        for nr=1:N_de_steps
            nn = Slice + NSlices * (nr-1);
            filename2= fullfile(pathname, sprintf(sprintf('%%s%%0%dd%%s', digitNo), fName(1:(idx - 1)), nn, ext));
            [tmp, infoMx]=dicom_read_singlefile(filename2, 1);
            DTI.dataAy(:,:,Slice ,nr)=tmp.dataAy;
            dE_direction = extract_DE_dirs(infoMx);
            DE_dirs(:,nr)= dE_direction;
        end
    end
else
    for nr = 1 : N_de_steps
        filename2 = create_full_filename(pathname,filename, nr);
        [DTI_volume, infoMx, msg]=dicom_read_singlefile(filename2, 1);    % this image will contain the whole volume
        for Slice= 1 : NSlices
            %%%% for mosaic images select the submatrix
            ColNumber= fix(Slice/NumberOfImagesInRaw - 0.00001) +1;
            RowNumber = Slice - (ColNumber-1)*NumberOfImagesInRaw;
            ColStart = 1+ sizeImOrig(1)*(ColNumber-1);
            RowStart = 1 + sizeImOrig(2)*(RowNumber-1);
            DTI.dataAy(:,:,Slice ,nr) = DTI_volume.dataAy(ColStart : ColStart+ sizeImOrig(1)-1 ,RowStart : RowStart + sizeImOrig(2) -1, : );
        end
        dE_direction = extract_DE_dirs(infoMx);
        DE_dirs(:,nr)= dE_direction;
    end
end

% update mrstruck info-fields to the proper format
DTI.dim3 ='size_z';
DTI.dim4 ='size_t';

% save structure
clear tmp DTI_volume
newFile = strcat(pathname,'DTI_',fName(1:idx(end)-3),'.mat');
mrstruct_write(DTI,newFile);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% start function  extract_DE_dirs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dE_direction = extract_DE_dirs(infoMx) %, ImageOrientationPatient)

for i=1:length(infoMx.Private_0029_1010),
    if(strcmp(infoMx.Private_0029_1010(i).name,'DiffusionGradientDirection' ) )
        iDiffusionGradientDirection =i;
        strucDiffusionGradientDirection= infoMx.Private_0029_1010(i);
    end
end

if isempty(strucDiffusionGradientDirection.item)
    dE_direction = [0 0 0];
else
    dE_direction = str2num([strucDiffusionGradientDirection.item(1).val strucDiffusionGradientDirection.item(2).val strucDiffusionGradientDirection.item(3).val]);
    %dE_direction= RotateDeGradients(tmp, ImageOrientationPatient);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of function extract_DE_dirs



%%%%%%%%% start function create_full_filename %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filename2 = create_full_filename(pathname,filename, filenumber )

if filenumber >=1000
    filename2 = [pathname, filename(1:length(filename)-8 ), num2str(filenumber),'.dcm'] ;
elseif filenumber >=100
    filename2 = [pathname, filename(1:length(filename)-8 ),'0', num2str(filenumber),'.dcm'] ;
elseif   filenumber >=10
    filename2 = [pathname, filename(1:length(filename)-8 ),'00', num2str(filenumber),'.dcm'] ;
else
    filename2 = [pathname, filename(1:length(filename)-8 ),'000', num2str(filenumber),'.dcm'] ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of function create_full_filename