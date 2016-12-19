% ReadDicomFilesDTIMZv5     for FR-Siemens sequences
% last update 030527 - interpolation to matrix 128x128
% ika 030227
% res = ReadDicomFilesDTI_MZv5(NSlices, N_de_steps, DTE_mode, NumberB0scans , tresholdAbs , b_factor, DO_interpolation, fileNumberingFlag)
% example:  res = ReadDicomFilesDTI_MZv5(4, 17, 'DE15', 2,  130 , 1000, 0, 'NonCont')
% based on ReadDicomCalcTensor4
% ultimative solution for GE and siemens data!
% ika 050413
% v.2 changed order of slice and DE loops
%
% v.3 ika050517 can work with huge files - does all slice per slice
% 050614 added option to roted the gradient directions properrly
%
%
% ika 050615    add Viewing tools   
%               add option for mosaic
% ika 051104: for motion and disctortion correction first image with comment
%             "reference" has to be removed
% ika 060405 - improved for mosaic images - now much faster!



% ! see %%%NB!!!!!%%% - has to be fixed - needs synchonisation with fiber TR




function [res, errStr] = ReadDicomFilesDTI_MZv5(NSlices, N_de_steps, DTE_mode, NumberB0scans , tresholdAbs , b_factor, DO_interpolation, fileNumberingFlag, vHd)


res= []; errStr= ''; % bwk
% errStr= strcat(mfilename, '(error): testError - has to be correctly implemented ')  % bwk
if ~exist('vHd') || ~ishandle(vHd)
    vHd= [];
end

errStr= 'calculation started';

DEBUG = 0; % set 1 for debug output

DoSingleSliceSAve = 'YES'; %'YES'; 'No '

[filename, pathname] = uigetfile('*.dcm', 'Select a file corresponding to the proper scan # ');

% go to 1st image
if isempty(vHd)
     disp(errStr);        
else
     set(vHd, 'string', errStr);
end
drawnow;
tic

filename1 = create_full_filename(pathname,filename, 1 ); % 1st image
filename2 = create_full_filename(pathname,filename, 2 ); % 2nd image

[tmp1, infoMx, msg]=dicom_read_singlefile(filename1, 1);

% check if the 1st image is a reference one
if isfield(infoMx,'ImageComments')
    if (~isempty(findstr(infoMx.ImageComments,'Reference')) )
        'Reference image is removed, Nr & N_bo reduced by 1'
        infoMx.ImageComments
        NumberB0scans =NumberB0scans-1;
        N_de_steps =N_de_steps-1;
        [tmp1, infoMx, msg]=dicom_read_singlefile(filename2, 1);
        deltaNumber = 1; % start from 2nd image
        ReferenceImageFlag ='Yes';
    else
        ReferenceImageFlag ='No';
        deltaNumber = 0; % increment for case if numeration  is uncontinuous as by stamps and MZ reco
    end
else
    ReferenceImageFlag ='No';
    deltaNumber = 0; % increment for case if numeration  is uncontinuous as by stamps and MZ reco
end

if ~exist('tresholdAbs') || isempty(tresholdAbs) % bwk
    tresholdAbs= 'def';                             % bwk
end                                             % bwk

if ~exist('fileNumberingFlag')
    fileNumberingFlag ='Cont'; %Continues file Numbering 
end;

for i=1:56, 
    if(strcmp(infoMx.Private_0029_1010(i).name,'NumberOfImagesInMosaic' ) ), i; strucNumberOfImagesInMosaic=infoMx.Private_0029_1010(i) ; end; 
end
wrongSlices = 0;
if(isempty(strucNumberOfImagesInMosaic.item)), 
    NumberOfImagesInMosaic = 1; 
    NumberOfImagesInRaw =1; % by mosaic is always NxN pictures 
else 
    NumberOfImagesInMosaic = str2num(strucNumberOfImagesInMosaic.item(1).val); 
    tmp = fix(sqrt(NumberOfImagesInMosaic));
    while (tmp^2 < NumberOfImagesInMosaic ), 
        tmp = tmp +1;
    end
    if ( tmp^2 >= NumberOfImagesInMosaic)
        NumberOfImagesInRaw =tmp;
    end
    if NumberOfImagesInMosaic < NSlices || NumberOfImagesInMosaic > NSlices
        NSlices = NumberOfImagesInMosaic; % slices given in dti_calulation gui does not agree with existing in mosaic images, correction!
        wrongSlices = 1;
    end
end

sizeImOrig =size(tmp1.dataAy); 
sizeImOrig(1) = sizeImOrig(1)/NumberOfImagesInRaw ;sizeImOrig(2)= sizeImOrig(2)/NumberOfImagesInRaw ; % size of each subscell in mosaic
sizeIm =sizeImOrig;

if(DO_interpolation >0)
    sizeIm =fix(sizeIm/2);
end


tmp1 =zeros( sizeIm(1),sizeIm(2) , N_de_steps);

if (strcmp( ReferenceImageFlag, 'Yes') )   
    if (NumberOfImagesInMosaic == 1), % noMosaic - old case 
        deltaNumber = NSlices; % scip 1st volume
    else
        deltaNumber = 1; % scip 1st mosaik image
    end
else
    deltaNumber = 0; % increment for case if numeration  is uncontinuous as by stamps and MZ reco
end

missed_filenumbers = []; % list off skiped file numbers

EigenVal1= [];
Max_error= [];
sq_Error= [];
eigVect1= [];
b0_image= [];


% ********************************************************************
% 
%     case of current (summer 2006) FR_siemens mosaic DTI data 
% 
% ********************************************************************

% if (NumberOfImagesInMosaic ~= 1)  save all data once in a float file

if (NumberOfImagesInMosaic ~= 1)
    index2 = strfind(filename, '_');
    fname =[pathname,'DWI_MZ_volumes_',filename(1:index2(end)-1), '.float' ];
    fid2=fopen(fname,'w+');
    if fid2<=2,
        disp('File could not be created');
    end
    for nr=1:N_de_steps
        outStr= sprintf('Prepare Data ... (%d/%d)', nr, N_de_steps);
        if isempty(vHd)
            disp(outStr);        
        else
            set(vHd, 'string', outStr);
            drawnow;
        end
        nn =  nr + deltaNumber;
        filename2 = create_full_filename(pathname,filename, nn );
        tmp_count = 0;
        while ( exist(filename2,  'file' )== 0)
            if (strcmp(fileNumberingFlag , 'Cont'))
                text =[' the next file  ', filename2, ' is not found, panic!']
                return
            end
            missed_filenumbers = [missed_filenumbers, nn];
            tmp_count =  tmp_count +1;
            nn = nn +1;
            filename2 = create_full_filename(pathname,filename, nn );
        end
        deltaNumber = deltaNumber + tmp_count;
        [DTI_volume, infoMx, msg]=dicom_read_singlefile(filename2, 1);    % this image will contain the whole volume
        
        for Slice= 1:NSlices
            %%%% for mosail images select the submatrix
            ColNumber= fix(Slice/NumberOfImagesInRaw - 0.00001) +1;
            RowNumber = Slice - (ColNumber-1)*NumberOfImagesInRaw;
            ColStart = 1+ sizeImOrig(1)*(ColNumber-1);
            RowStart = 1 + sizeImOrig(2)*(RowNumber-1);
            
            if(DEBUG)
                tmpX2 = DTI.dataAy;    
                tmpX2(ColStart : ColStart+ sizeImOrig(1)-1 ,RowStart : RowStart + sizeImOrig(2) -1, : ) = 0*tmpX2(ColStart : ColStart+ sizeImOrig(1)-1 ,RowStart : RowStart + sizeImOrig(2) -1, : );
                figure(Slice)
                imagesc(tmpX2)
            end
            
            volume(:,:,Slice) = DTI_volume.dataAy(ColStart : ColStart+ sizeImOrig(1)-1 ,RowStart : RowStart + sizeImOrig(2) -1, : );
            %%%%
            
        end   %Slice= 1:NSlices
        
        if(DO_interpolation >0)
            %         sizeIm =size(DTI.dataAy);
            tmp = interpft(volume , sizeIm(1) , 1); 
            volume  = interpft(tmp ,sizeIm(2) , 2); 
        end
        
%         if(DEBUG)
%             infoMx.ImageComments
%         end
        
        time = infoMx.AcquisitionTime;
        if nr == 1
            hour1 = str2double(time(1,1:2))*60*60;
        end
        hours = (str2double(time(1,1:2))*60*60-hour1); % in seconds % if in milliseconds: * 1000
        minutes = str2double(time(1,3:4))*60; % in seconds % if in milliseconds: * 1000
        seconds = str2double(time(1,5:6)); % in seconds % if in milliseconds: * 1000
        micro = str2double(time(1,8:end))/1000000; % in seconds % if in milliseconds: /1000
        AcquisitionTime(1,nr) =  hours + minutes + seconds + micro; % in seconds, each image has his own time!
        SliceLocation(1,nr) =  (infoMx.SliceLocation);
        ImagePositionPatient =(infoMx.ImagePositionPatient); % [3x1 double]
        ImageOrientationPatient =(infoMx.ImageOrientationPatient); % [6x1 double]
        if isfield(infoMx,'ImageComments')
            ImageComments =infoMx.ImageComments; % contains gradients orienbtation
        else
            ImageComments = '';
        end

        % have to sort  acording slice position
        % !!! now can not sort - if needed has to be done later
        [ax1, ay1] =sort(SliceLocation(:,1));
        %     %tmp1(:,:,ay1 ,:)=tmp1(:,:,1:end,:);
        %     tmp1=tmp1(:,:,ay1,:);
        SliceTimeOrder=ay1;
        
%        [dE_direction, AcquisitionMatrixText,  b_value(nr)] = extract_DWI_parameters(infoMx, ImagePositionPatient);
        [dE_direction, AcquisitionMatrixText,  b_value(nr)] = extract_DWI_parameters(infoMx, ImageOrientationPatient);
        DE_direction(1,nr,:) = dE_direction;
        
        fseek(fid2, 0, 'eof');
        resWrite = fwrite(fid2,volume,'float32');
    end % for nr=1:N_de_steps
    
    text =['save(''', pathname,'DWI_infoMZ_mosaic_',filename(1:index2(end)-1), '.mat'',  ''DE_direction'',  ''AcquisitionTime'', ''SliceLocation'');' ]; % SliceTimeOrder' ]
    eval(text) 
    fclose(fid2); % close write file
    
    % now test - read and save per slices
    DEs = squeeze(DE_direction);
    nulls_ind = find(DEs(:,1) == 0 & DEs(:,2) == 0 & DEs(:,3) == 0);
    nulls_length = length(nulls_ind);
    DTI_slice = DTI_volume;
    DTI_slice.dim3= 'size_t';
    
        %%%%%%%%%%%%%%%%%%%%% READ DTA FROM FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    index2 = strfind(filename, '_');
    fname =[pathname,'DWI_MZ_volumes_',filename(1:index2(end)-1), '.float' ];
    fid3 = fopen(fname,'r');  % open data file for read
    offs=0;
    for Slice= 1:NSlices
        offs = (Slice-1)*4*sizeIm(1)*sizeIm(2);
        for nr=1:N_de_steps
            fseek(fid3,offs,-1);
            a=fread(fid3,[sizeIm(1),sizeIm(2) ],'float32');  % read  data from file
            DataSlice(:,:,nr) =a;
            offs = offs+4*sizeIm(1)*sizeIm(2)*NSlices;
        end % for nr=1:N_de_steps
        
            %%%%%%%%%%%%%%%%%%%%% DONE READ DTA FROM FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        DTI_slice.dataAy = DataSlice;
        b0s = DTI_slice.dataAy(:,:,nulls_ind);
        DTI_slice.dataAy(:,:,nulls_ind) = [];
        DTI_slice.dataAy = cat(3,b0s,DTI_slice.dataAy);
        %DWI =DTI_slice.dataAy;
        if(strcmp(DoSingleSliceSAve , 'YES'))
            eval(['save(''DWI_', filename(1:index2(end)-1),'_', num2str(Slice), ''', ''DTI_slice'');' ]);     
        end
        [D_trace, FA, RA, EigenVal, eigVect, M_error, sqErr] = do_dti_calculation(DE_direction, b_factor, DataSlice, sizeIm , tresholdAbs )  ;
        
        if isempty(EigenVal1)
            siz= size(EigenVal);
            EigenVal1= zeros(siz(1), siz(2), NSlices, 3);
            eigVect1= zeros(siz(1), siz(2), NSlices, 3, 3);
            Max_error= zeros(siz(1), siz(2), NSlices);
            sq_Error= zeros(siz(1), siz(2), NSlices);
            b0_image= zeros(siz(1), siz(2), NSlices);
        end
    
        EigenVal1(:,:,Slice,:) = EigenVal; %ika2k1206
        eigVect1(:,:,Slice,:,:) = eigVect;
        b0_image(:,:,Slice) = DataSlice(:,:,1);

        Max_error(:,:,Slice) = M_error;
        sq_Error(:,:,Slice)= sqErr;
        res = toc;
        outputS=sprintf('done for slice %d/%d (time used: %dmin %ss)', Slice, NSlices, floor(res/60), num2str(round((res - 60*floor(res/60))*10)/10));
        if isempty(vHd)
            disp(outputS);        
        else
            set(vHd, 'string', outputS);
            drawnow;
        end
    end % for Slice= 1:NSlices
    fclose(fid3); % close read from file
        
    b0_image_struc= mrstruct_init('volume', b0_image, DTI_volume);
    error_struc= mrstruct_init('volume', Max_error, b0_image_struc);
    eigenVal_struc = mrstruct_init('series3D', EigenVal1, b0_image_struc);
    eigenVec_struc = mrstruct_init('series3DEchos', eigVect1, b0_image_struc);

    FileName =[pathname,filename(1:index2(end)-1)];
    text = ['save(''', FileName,'_DTD.mat'', ''b0_image_struc'', ''eigenVal_struc'', ''eigenVec_struc'', ''error_struc'');'];
    eval(text)
    
    % done save dti structures
    if wrongSlices == 1
        outStr = 'Calculation finished. Wrong no of slices was given! Check output!!';
        'Calculation finished. Wrong no of slices was given! Check output!!', NSlices
    else
        outStr = 'calculation finished';
    end
    if isempty(vHd)
        disp(outStr);
    else
        set(vHd, 'string', outStr);
        drawnow;
    end
    res =toc;
    return
    
end %if (NumberOfImagesInMosaic ~= 1)

% ********************************************************************
% 
%     FINISHED case of current (summer 2006) FR_siemens mosaic DTI data 
% 
% ********************************************************************


	% 	
	% 	% ********************************************************************
	% 	% 
	% 	%     below is the code for the old FR_siemens DTI data , actually only
	% 	%     single_slice data will go through
	% 	% 
	% 	% ********************************************************************
	% 	
	% 	
	% 	% open file to write in "old Format" - i.e. for each slice all DE steps,
	% 	% than next slice .. etc.
	% 	index2 = findstr(filename, '_');
	% 	fname =[pathname,'DWI_MZ_',filename(1:index2(end)-1), '.float' ];
	% 	fid=fopen(fname,'a+');
	% 	if fid<=2,
	%         disp('File could not be created');
	% 	else
	%         ;
	% 	end
	% 	
	% 	
	% 	for Slice= 1:NSlices
	%         for nr=1:N_de_steps
	%             if (NumberOfImagesInMosaic == 1), % noMosaic - old case
	%                 nn =  nr + (Slice-1)*N_de_steps + deltaNumber;
	%                 filename2 = create_full_filename(pathname,filename, nn );
	%                 tmp_count = 0;%
	%                 while ( exist(filename2,  'file' )== 0)
	%                     if (strcmp(fileNumberingFlag , 'Cont') )
	%                         text =[' the next file  ', filename2, ' is not found, panic!']
	%                         return
	%                     end
	%                     missed_filenumbers = [missed_filenumbers, nn];
	%                     tmp_count =  tmp_count +1;
	%                     nn = nn +1;
	%                     filename2 = create_full_filename(pathname,filename, nn );
	%                 end
	%                 deltaNumber = deltaNumber + tmp_count ;
	%                 [DTI, infoMx, msg]=dicom_read_singlefile(filename2, 1);
	% 	
	%             else    % Mosaic - new case
	%                 nn =  nr +  deltaNumber; % I accume that the file numeration is correct and continious now!
	%                 filename2 = create_full_filename(pathname,filename, nn );
	%                 [DTI, infoMx, msg]=dicom_read_singlefile(filename2, 1);
	%                 %%%% for mosaic images select the submatrix
	%                 ColNumber= fix(Slice/NumberOfImagesInRaw - 0.00001) +1;
	%                 RowNumber = Slice - (ColNumber-1)*NumberOfImagesInRaw;
	%                 ColStart = 1+ sizeImOrig(1)*(ColNumber-1);
	%                 RowStart = 1 + sizeImOrig(2)*(RowNumber-1);
	% 	
	%                 if(DEBUG)
	%                     tmpX2 = DTI.dataAy;
	%                     tmpX2(ColStart : ColStart+ sizeImOrig(1)-1 ,RowStart : RowStart + sizeImOrig(2) -1, : ) = 0*tmpX2(ColStart : ColStart+ sizeImOrig(1)-1 ,RowStart : RowStart + sizeImOrig(2) -1, : );
	%                     figure(Slice)
	%                     imagesc(tmpX2)
	%                 end
	% 	
	%                 DTI.dataAy = DTI.dataAy(ColStart : ColStart+ sizeImOrig(1)-1 ,RowStart : RowStart + sizeImOrig(2) -1, : );
	%                 %%%%
	%             end
	% 	
	% 	%         if(DEBUG)
	% 	%             infoMx.ImageComments
	% 	%         end
	%             % &&&&&&&&&&&&&&&&&
	%             if(DO_interpolation >0)
	%                 %         sizeIm =size(DTI.dataAy);
	%                 tmp = interpft(DTI.dataAy , sizeIm(1) , 1);
	%                 DTI.dataAy = interpft(tmp ,sizeIm(2) , 2);
	%             else
	%                 ;
	%             end
	%             %%%%%%%%%%%        save importand parameters
	%             time = infoMx.AcquisitionTime;
	%             if nr == 1
	%                 hour1 = str2double(time(1,1:2))*60*60;
	%             end
	%             hours = (str2double(time(1,1:2))*60*60-hour1); % in seconds % if in milliseconds: * 1000
	%             minutes = str2double(time(1,3:4))*60; % in seconds % if in milliseconds: * 1000
	%             seconds = str2double(time(1,5:6)); % in seconds % if in milliseconds: * 1000
	%             micro = str2double(time(1,8:end))/1000000; % in seconds % if in milliseconds: /1000
	%             AcquisitionTime(Slice,nr) =  hours + minutes + seconds + micro; % each image has his own time!
	%             SliceLocation(Slice,nr) =  (infoMx.SliceLocation);
	%             %%%                how to handle this????
	%             ImagePositionPatient =(infoMx.ImagePositionPatient); % [3x1 double]
	%             ImageOrientationPatient =(infoMx.ImageOrientationPatient); % [6x1 double]
	%             if isfield(infoMx,'ImageComments')
	%                 ImageComments =infoMx.ImageComments; % contains gradients orienbtation
	%             else
	%                 ImageComments = '';
	%             end
	%             %        [dE_direction, AcquisitionMatrixText,  b_value(nr)] = extract_DWI_parameters(infoMx, ImagePositionPatient);
	%             [dE_direction, AcquisitionMatrixText,  b_value(nr)] = extract_DWI_parameters(infoMx, ImageOrientationPatient);
	%             DE_direction(Slice,nr,:) = dE_direction;
	% 	
	%             tmp1(:,:,nr)=DTI.dataAy; %from v.3
	%         end  % for nr=1:N_de_steps
	% 	
	%         % save raw data for the slice in float format
	%         %     index2 = findstr(filename, '_');  %ika070723  what was sence of these 2 lines?
	%         %     fname =[pathname,'DWI_MZ_',filename(1:index2(end)-1), '.float' ];%ika070723 ...these 2 lines??
	%         resWrite = fwrite(fid,tmp1,'float32');
	%         % done save for the current slice
	%         if(strcmp(DoSingleSliceSAve , 'YES'))
	%             DWI = DTI;
	%             DWI.dataAy = tmp1;
	%             DWI.dim3 = 'size_t';
	%             text = ['save(''DWI_', filename(1:index2(end)-1),'_', num2str(Slice), ''', ''DWI'');' ];
	%             eval(text)
	%         end
	% 	
	%         % have to sort  acording slice position
	%         % !!! now can not sort - if needed has to be done later
	%         [ax1, ay1] =sort(SliceLocation(:,1));
	%         %     %tmp1(:,:,ay1 ,:)=tmp1(:,:,1:end,:);
	%         %     tmp1=tmp1(:,:,ay1,:);
	% 	
	%         SliceTimeOrder=ay1;
	%         % done sort
	% 	
	%         DTI.dataAy=tmp1;
	%         DTI.dim3 ='size_z'; % set mrstruct properly
	%         DTI.dim4 ='size_t'; % set mrstruct properly
	%         if(DO_interpolation >0)
	%             DTI.vox(1)= 2*DTI.vox(1);
	%             DTI.vox(2)= 2*DTI.vox(2);
	%         end
	%         toc
	% 	
	%         [D_trace, FA, RA, EigenVal, eigVect, M_error, sqErr] = do_dti_calculation(DE_direction, b_factor, tmp1,  sizeIm , tresholdAbs );
	% 	
	%         if isempty(EigenVal1)
	%             siz= size(FA);
	%             EigenVal1= zeros(siz(1), siz(2), NSlices, 3);
	%             eigVect1= zeros(siz(1), siz(2), NSlices, 3, 3);
	%             Max_error= zeros(siz(1), siz(2), NSlices);
	%             sq_Error= zeros(siz(1), siz(2), NSlices);
	%             b0_image= zeros(siz(1), siz(2), NSlices);
	%         end
	% 	
	%         EigenVal1(:,:,Slice,:)= EigenVal; %ika2k1206
	%         Max_error(:,:,Slice) = M_error;
	%         sq_Error(:,:,Slice)= sqErr;
	%         eigVect1(:,: ,Slice, :,:) =eigVect;
	%         b0_image(:,:,Slice)= tmp1(:,:,1 );
	% 	
	%         outputS=sprintf('done for slice %d/%d (time used: %dmin %ss)', Slice, NSlices, floor(res/60), num2str(round((res - 60*floor(res/60))*10)/10));
	%         if isempty(vHd)
	%             disp(outputS);
	%         else
	%             set(vHd, 'string', outputS);
	%             drawnow;
	%         end
	%         % % %     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%         % % %     %
	%         % % %     %       DONE DTI calculation
	%         % % %     %
	%         % % %     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%         Slice
	% 	end; % Slice
	% 	fclose(fid); % close write file
	% 	text =['save(''', pathname,'DWI_infoMZ2_',filename(1:index2(end)-1), '.mat'', ''DE_direction'', ''AcquisitionTime'', ''SliceLocation'');' ] % SliceTimeOrder' ]
	% 	eval(text)
	% 	
	% 	% b0_image_struc= mrstruct_init('volume', DTI_volume, b0_image);  %%
	% 	% ika061011   - bjoern I don't understand that
	% 	%
	% 	b0_image_struc= mrstruct_init('volume',b0_image, DTI );
	% 	error_struc= mrstruct_init('volume', Max_error, b0_image_struc);
	% 	eigenVal_struc = mrstruct_init('series3D', EigenVal1, b0_image_struc);
	% 	eigenVec_struc = mrstruct_init('series3DEchos', eigVect1, b0_image_struc);
	% 	
	% 	FileName =[pathname,filename(1:index2(end)-1)];
	% 	text = ['save ', FileName,'_DTD.mat',  ' b0_image_struc eigenVal_struc eigenVec_struc error_struc'];
	% 	eval(text)
	% 	
	% 	% done save dti structures
	% 	
	% 	outStr = 'calculation finished';
	% 	if isempty(vHd)
	%         disp(outStr);
	% 	else
	%         set(vHd, 'string', outStr);
	%         drawnow;
	% 	end
	% 	
	% 	res =toc
	% 	return

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   END OF MAIN
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% local functions
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of function create_full_filename
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% function  extract_DWI_parameters

%function [dE_direction, AcquisitionMatrixText,  b_value] = extract_DWI_parameters(infoMx, ImagePositionPatient)
function [dE_direction, AcquisitionMatrixText,  b_value] = extract_DWI_parameters(infoMx, ImageOrientationPatient)


for i=1:length(infoMx.Private_0029_1010),
    if(strcmp(infoMx.Private_0029_1010(i).name,'DiffusionGradientDirection' ) )
        iDiffusionGradientDirection =i;
        strucDiffusionGradientDirection= infoMx.Private_0029_1010(i);
    end
    if(strcmp(infoMx.Private_0029_1010(i).name,'DiffusionDirectionality' ) ),
        iDirectionality=i; strucDirectionality=infoMx.Private_0029_1010(i);
    end
    if(strcmp(infoMx.Private_0029_1010(i).name,'B_value' ) ),
        iB_value=i; strucB_value=infoMx.Private_0029_1010(i);
    end
    %   if(strcmp(infoMx.Private_0029_1010(i).name,'DiffusionDirectionality' ) ), i, strucDiffusionDirectionality =infoMx.Private_0029_1010(i) , end;
    if(strcmp(infoMx.Private_0029_1010(i).name,'NumberOfImagesInMosaic' ) )
        strucNumberOfImagesInMosaic=infoMx.Private_0029_1010(i);
    end
    if(strcmp(infoMx.Private_0029_1010(i).name,'AcquisitionMatrixText' ) )
        strucAcquisitionMatrixText =infoMx.Private_0029_1010(i);
    end
    %   if(strcmp(infoMx.Private_0029_1010(i).name,'SliceNormalVector' ) ), i, strucSliceNormalVector=infoMx.Private_0029_1010(i), end;
    %   if(strcmp(infoMx.Private_0029_1010(i).name,'TimeAfterStart' ) ), i, strucTimeAfterStart=infoMx.Private_0029_1010(i) ,end;
end

%matrix size as text
AcquisitionMatrixText = strucAcquisitionMatrixText.item(1);

%         %NumberOfImagesInMosaic - can be empty
%         if(isempty(strucNumberOfImagesInMosaic.item)), NumberOfImagesInMosaic = 1; else NumberOfImagesInMosaic = str2num(strucNumberOfImagesInMosaic.item(1).val); end
%
% diffusion specific parameters
if(isempty(strucB_value.item)), b_value = 99999999999999999999, else b_value = str2num(strucB_value.item(1).val); end %if there now b-val in the corresponding field give just come strange value

%DE direction
if isempty(strucDiffusionGradientDirection.item) && strcmp(deblank(strucDirectionality.item(1,1).val), 'DIRECTIONAL')
    a =strfind(ImageComments,'b=');
    if(~isempty(strfind(ImageComments,'b=0,')))
        dE_direction(:) = [0 0 0]; % b0 scan
    elseif (~isempty(a))
        b =strfind(ImageComments,':');
        tmp = str2num(ImageComments(b+1:end));
        dE_direction(:) = RotateDeGradients(tmp, ImageOrientationPatient);
    else
        a2 =strfind(ImageComments,'DW: ');
        if(~isempty(a2))
            % DE_direction(Slice,nr,:) = str2num(ImageComments(a+4:end));
            tmp = str2num(ImageComments(a2+4:end));
            dE_direction(:) = rotateDeGradients(tmp, ImageOrientationPatient);
        elseif (~isempty(strfind(ImageComments,'No DW')) )
            dE_direction = [0 0 0]; % b0 scan
        else
            dE_direction = [10 10 10] % set norm >1 to mark non commented scan
        end
    end
elseif isempty(strucDiffusionGradientDirection.item) && strcmp(deblank(strucDirectionality.item(1,1).val), 'NONE')
    dE_direction = [0 0 0];
else
%    tmp = str2num([strucDiffusionGradientDirection.item(1).val strucDiffusionGradientDirection.item(2).val strucDiffusionGradientDirection.item(3).val]);
%    dE_direction= RotateDeGradients(tmp, ImageOrientationPatient);
   dE_direction = str2num([strucDiffusionGradientDirection.item(1).val strucDiffusionGradientDirection.item(2).val strucDiffusionGradientDirection.item(3).val]);
end

%%%NB!!!!!%%%
dE_direction(2) = -dE_direction(2);% 31.10.03 - why I have done that????
dE_direction(1) = -dE_direction(1);% 31.10.03 - why I have done that????
%%%NB!!!!!%%%

%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%% function  do_dti_calculation
%
%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D_trace, FA, RA, EigenVal, eigVect, M_error, sqErr] = do_dti_calculation(DE_direction, b_factor, singleSliceDWI_set, sizeIm, tresholdAbs )

D_trace= []; FA= []; RA= [];  EigenVal= []; eigVect= []; M_error= []; sqErr= [];
DTE_mode  = DE_direction ; % transfer DE grad directions read from dicom comment field
%% now time to created the b-matrix
DE = squeeze(DTE_mode(1,:,:)); % !! NB now use the data for the 1st slice, if motion correction is on has to be fixed ika 050413
B_tensor =zeros(3,3,length(DE));
DE_size = size(DE);
for i = 2:DE_size(1)
    tmp=DE(i,:)'*DE(i,:);
    if sum(sum(abs(tmp) )) >0, % ika 041118
        B_tensor(:,:,i)= (b_factor/trace(tmp))*tmp;end
end
nsteps =length(B_tensor);
X1 = [  reshape(B_tensor(1,1,:), 1, nsteps , 1)]; %1st diag element
X2 = [  reshape(B_tensor(2,2,:), 1, nsteps , 1)]; %2nd diag element
X3 = [  reshape(B_tensor(3,3,:), 1, nsteps , 1)]; %3rd diag element
X4 = [  reshape(B_tensor(1,2,:), 1, nsteps, 1)];
X5 = [  reshape(B_tensor(1,3,:), 1, nsteps, 1)];
X6 = [  reshape(B_tensor(2,3,:), 1, nsteps, 1)];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   do dti calculation for the paticular slice
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mask =zeros(sizeIm(1), sizeIm(2));
treshold = mean(mean(singleSliceDWI_set(:,:,1 )));  % NB now maks will depend from the slice - for low slices that can be bad
D_trace= mask;
D_iso= mask;

for i=1:sizeIm(1) % vertical
    for j=1:sizeIm(2) % horizontal
        if singleSliceDWI_set(i,j,1) > tresholdAbs %treshold_factor *treshold
            mask(i,j)=1;
        end
    end
end

% reserviere speicher
D_trace= zeros(sizeIm(1:2));
eigVect=zeros(sizeIm(1), sizeIm(2), 3,3);
difComp=zeros(sizeIm(1), sizeIm(2), 3,3);
M_error= zeros(sizeIm(1:2));
sqErr= zeros(sizeIm(1:2));


Y1 =singleSliceDWI_set;
for x_coor =1: sizeIm(2)
    for y_coor =1: sizeIm(1)
        Y  = [];	% clear
        if (mask(y_coor, x_coor)) ~= 0;
            tempY =(Y1(  y_coor, x_coor, :));
            itemp =find(tempY <= 0);
            if(length(itemp) >=1),
                %     'warning! signal is zero or negative',tempY(itemp),  tempY(itemp)= 20.0222; x_coor, y_coor,%  itemp,
                tempY(itemp)= 20.0222;
            end
            y=squeeze( sum( sum( log(tempY ), 2 ),1 ) );  % 2:5 - from 2 to 5th image
            Y= [Y y'];
            anisoX=[ones(size(X1')), X1', X2', X3',2*X4', 2*X5', 2*X6'];
            aniso_a=anisoX\(Y');
            YY = anisoX*aniso_a;
            MaxErr = max(abs(YY- Y'));

            d_tens= [aniso_a(2), aniso_a(5) ,aniso_a(6) ;aniso_a(5) ,aniso_a(3) ,aniso_a(7); aniso_a(6) ,aniso_a(7) ,aniso_a(4) ] ; % check off_diaf elem order!!!!
            [vv, dd] = eig(-d_tens) ;
            %'vv -eigen-vector, dd - eigen-value';
            D_trace(y_coor, x_coor ) = sum(sum(dd)) * 1/3; % ika 2k0121 added 1/3 factor
            eigVect(y_coor, x_coor , :,:) =vv(:,:);
            difComp(y_coor, x_coor , :,:) =dd(:,:);
            % D_raw(y_coor, x_coor,:,: ) = -d_tens;  % ika9041130
            M_error(y_coor, x_coor ) = MaxErr;
            sqErr(y_coor, x_coor ) = mean((YY- Y').^2);
        end %if
    end % for
end % for

EigenVal= sum(difComp, 3);
dd_1= EigenVal(:,:,1).^2 + EigenVal(:,:,2).^2 + EigenVal(:,:,3).^2;
l_mean= (EigenVal(:,:,1) + EigenVal(:,:,2) + EigenVal(:,:,3))/3;
dd_3= (EigenVal(:,:,1)-l_mean).^2 +(EigenVal(:,:,2)-l_mean).^2 +(EigenVal(:,:,3)-l_mean).^2 ;
i6= find(dd_3 ~= 0);
i7= find(l_mean ~= 0);
FA = zeros(size(dd_3));
RA = FA;
FA(i6)=sqrt(3/2)* sqrt(dd_3(i6)./dd_1(i6));
RA(i7)=sqrt(1/3)* sqrt(dd_3(i7))./l_mean(i7);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%       DONE DTI calculation
%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% end of file %%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
