% interface program for batch program
% roi log
% Author: Susanne Schnell
% Windows 22.02.2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = dti_logroi_ui(P)

logFName= fullfile(P.logpath{1}, P.newfilename);
dtd = dtdstruct_read(P.dtdname{1});
comStr= dtdstruct_query(dtd, 'getComStr');
ROI = maskstruct_read(P.roiname{1});
ROInames = maskstruct_query(ROI,'maskNames');


if isfield(P.status,'ValsRoi')
    number = P.status.ValsRoi.mask1;
    if number(2) == Inf
        number = 1 : 1 : numel(ROInames);
    end
    if isempty(dtd.b0_image_struc.patient)
        patientStr= 'dummyName';
    else
        patientStr= str2filename(dtd.b0_image_struc.patient);
    end
    % hole alle Bilddaten und eigenwerte  
    for i= 1:length(comStr)
        dataVal= dtdstruct_query(dtd, strcat('get', comStr{i}));
        paramCell{i}= dataVal.dataAy;
    end
    
    % Fuer alle nicht leeren ROIS 
    for i = 1:length(number)
        roiName= ROInames{number(i)};
        [idx, errStr]= maskstruct_query(ROI, 'getMaskIdx', roiName);
        if ~isempty(errStr)
            disp(errStr);
            return
        end
        if ~isempty(idx)
            logFName= fullfile(P.logpath{1}, str2filename(sprintf(P.newfilename, patientStr, roiName), 'extFileName'));
            file_ID= fopen(logFName, 'wt');
            if file_ID < 0
                logFName= fullfile(P.logpath, str2filename(sprintf(P.newfilename, patientStr, roiName)));
                file_ID= fopen(logFName, 'wt');
                if file_ID < 0
                    disp(strcat('Could not create file ''', logFName, ''''));
                    return
                end
            end
            dataMy= zeros(length(idx), length(comStr));
            for j= 1:length(comStr)
                fprintf(file_ID, ' %s ', comStr{j});
                dataMy(:, j)= reshape(paramCell{j}(idx), [length(idx) 1]);
            end
            valStrMy= reshape(sprintf('%20.13g', dataMy')', [size(dataMy, 2)*20, size(dataMy, 1)])';
            %                num2str(dataMy, 10);
            valStrMy(:, end)= double(sprintf('\n'));
            fprintf(file_ID, '\n%s\n', valStrMy');
            fclose(file_ID); %close file
        end
    end
    
else %log stats

    modNo= length(comStr);
    headLineStr= sprintf('Name\t DateOfBirth\t DataFileName\t ExaminationDate\t Roi_Name\t date_of_ROI_saving\t NumberOfPixels');
    
    for i= 1:modNo
        headLineStr= sprintf('%s\t mean(%s)\t std(%s)', headLineStr, comStr{i}, comStr{i});
    end
    
    if exist(logFName, 'file') == 2
        file_ID= fopen(logFName, 'r+t');
        if isempty(file_ID) || (file_ID < 0)
            disp(strcat('Could not reopen file ''', logFName, ''''));
            return
        end
        fseek(file_ID, 0, 'bof');
        lStr= fgetl(file_ID);
        fclose(file_ID);
        if ~strcmp(lStr, headLineStr)
            butStr= questdlg('The existing log file seems not to be compatible with the current dtdStrucht', ...
                'ROI logging', 'continue', 'truncate', 'abort','abort');
            if strcmp(butStr, 'continue')
                file_ID= fopen(logFName, 'at');
            elseif strcmp(butStr, 'truncate')
                file_ID= fopen(logFName, 'wt');
                fprintf(file_ID, '%s\n', headLineStr);
            else
                return;
            end
        else
            file_ID= fopen(logFName, 'at');
        end
    else
        file_ID= fopen(logFName, 'wt');
        if isempty(file_ID) || (file_ID < 0)
            disp(strcat('Could not create file ''', logFName, ''''));
            return
        end
        fprintf(file_ID, '%s\n', headLineStr);
    end
    
    if isempty(P.dtdname)
        fName= '<unknown>';
    else
        [p, fName, e]= fileparts(P.dtdname{1});
    end
    
    patientStr= dtdstruct_query(dtd, 'patient');
    if isempty(patientStr)
        patientStr= '<unknown>';
    end
    
    maskNames= maskstruct_query(ROI, 'maskNames');
    for roiI= 1:numel(maskNames)
        disp(['Logging ROI ',maskNames{roiI},' in ', logFName]);

        roiIdx= maskstruct_query(ROI, 'getMaskIdx', maskNames{roiI});
        fprintf(file_ID, '%s\t %s\t %s\t %s\t %s\t %s\t %d',  ...
            patientStr, '<unknown>', fName, '<unknown>', maskNames{roiI}, date, length(roiIdx));
        
        for i= 1:modNo
            tmpData = dtdstruct_query(dtd, strcat('get', comStr{i}));
            if isempty(roiIdx)
                fprintf(file_ID, '\t %s\t %s', ...
                    '<undef>', ...
                    '<undef>');
            elseif length(roiIdx) == 1
                fprintf(file_ID, '\t %s\t %s', ...
                    norm_numbers(tmpData.dataAy(roiIdx), 1, 8, 'E', 3), ...
                    '<undef>');
            else
                fprintf(file_ID, '\t %s\t %s', ...
                    norm_numbers(mean(tmpData.dataAy(roiIdx)), 1, 8, 'E', 3), ...
                    norm_numbers(std(tmpData.dataAy(roiIdx)), 1, 8, 'E', 3));
            end
        end
        fprintf(file_ID, '\n');
    end
    fclose(file_ID);
    disp(strcat('logROI: Map statistics are saved in ''', logFName, ''''));
end

out.files{1} = P.newfilename;