% read_dti_dicom: Read all dicom DTI data indepently from manufacturer
%
%       run data through 'dicom_copy_tool' first!!
%
%       [filename,status,err] = read_dti_dicom(path,list,flag_averages,fName,count_averages,vHd)
%
%       input: 
%       -'path' = path where dicom files are located
%       -'list' = list of filenames as structure
%       -'flag_averages' = if flag is 0 then no averages, else averages and fName has to
%           be given as well!
%       -'fName' = name of temporary file without endings ('_info.mat','_raw.bin')
%       -'count_averages' = actual number of average volume
%       -'vHd' = GUI handle
%       -> works also without input!
%
%       output:
%       -'filename' = filepath of saved data in ('_info.mat','_raw.bin')
%           files
%       -'status' = vector with length of 3 giving status if diffusion
%       information is complete: status(1) == 0 => no DEscheme, status(2) ==
%       0 => number of b0 images unknown, status(3) == 0 => b-value is
%       unknown; if one of these scalars equals one this specific
%       information is known and saved in '_info.mat'
%       -'err' = error string
%
%       Example usage:
%       path = 'E:/data_sets';
%       list = dir(strcat(path,'\*.dcm')); % here for PC
%       [filename,status,err] = read_dti_dicom(path,list,0,[])
%       --> no averages, no gui handle
%
%       function uses tho following functions:
%           dicom_read_singlefile
%           dw_data_admin
%           get_volume_geo
%           hm_create
%           hm_analyse
%           ...
%       
%       Please take care:
%       - not thoroughly tested yet!
%
% code and knowledge partly copied from fromer group member Kamil Il'yasov (University Kazan)
%__________________________________________________________________________
%		Susanne Schnell
%       Department of Diagnostic Radiology, Medical Physics
%       University Hospital Freiburg
%       Hugstetter Strasse 55
%       79106 Freiburg, Germany
%		08-12/2007
%
%-------------------- development -------------------------------------
%   AUTHORS:
%           Susanne Schnell (SuS)
%           Kamil Il'yasov

function [varargout] = read_dti_dicom(path, list, flag_averages,fName,count_averages, vHd)

if nargin > 6
    varargout{2} = 0;
    varargout{3} = 'Error in read_DTI_dicom: too many input arguments.';
    return;
end

if nargin == 0
    vHd = [];
    flag_averages = 0;
    fName = [];
    count_averages = [];
    path = uigetdir(pwd,'Select the directory containing the dicom files!'); % Matlab 2007b
    %path = 'E:\MR\DTI\DataSets\ThimoGrotz_15012008\DTI_slices69_vox2x2x2_matr104x104_TR11000_TE94_DEdirs61_8b0s'; % Matlab R14
    list = dir(fullfile(path,'*.dcm'));
    if isequal(list , 0) || isequal(path , 0)
        disp('User pressed cancel');
        return;
    else
        fprintf('User selected dcm files from %s', path);
    end
elseif nargin == 1
    vHd = [];
    flag_averages = 0;
    fName = [];
    count_averages = [];
    list = dir(fullfile(path,'*.dcm'));
    while isempty(list(1,1).name)
        disp('no dicoms in folder, please try again');
        path = uigetdir(pwd,'Select the directory containing the dicom files!'); % Matlab 2007b
            %path = uigetdir('','Select the folder containing the dicom files'); % Matlab R14
        list = dir(fullfile(path,'*.dcm'));
    end
elseif nargin == 2
    vHd = [];
    flag_averages = 0;
    fName = [];
    count_averages = [];
elseif nargin == 3
    count_averages = [];
    vHd = [];
    if flag_averages == 0
        fName = [];
    else
        varargout{2} = 0;
        varargout{3} = 'Error in read_DTI_Dicom: name of temporary data file missing!';
        return;
    end
elseif nargin == 4
    count_averages = [];
    vHd = [];
elseif nargin == 5
    vHd = [];
end

if ~ishandle(vHd)
    vHd= [];
end
if flag_averages == 0
    mes= 'data preperation started';
else
    mes = 'averaging started';
end
if isempty(vHd)
    disp(mes);
else
    set(vHd, 'string', mes);
    drawnow;
end

% write path in front of file names
for m = 1 : size(list,1)
    list(m,:).name = fullfile(path,list(m,:).name);
end

% function to read header information and put it in info file, only when
% this is the first of the averages or no averages at all
if flag_averages == 0
    try
        [information,user,fName1] = get_header_info(list,vHd);
        if isempty(fName)
            fName = fName1;
        end
        if isempty(vHd)
            disp('saving info file ...');
        else
            set(vHd,'string','saving info file ...');
            drawnow
        end
    catch
        mes = lasterr; 
        mes = strcat('Error in read_DTI_Dicom: ',mes);
        if isempty(vHd)
            disp(mes);
        else
            set(vHd,'string',mes);
        end
        varargout{2} = 0;
        varargout{3} = mes;
        return;
    end
    %create data and info file and save info data using dw_data_admin
    [res, errStr]= dw_data_admin('create',user.matrix,user.SliceNo,user.acquNo,fName);
    if res == 0
        if isempty(vHd)
            disp('Error in read_dti_dicom: files could not been created!');
            disp(errStr);
        else
            set(vHd,'string','Error in read_dti_dicom: files could not been created!');
            disp(errStr);
        end
        varargout{2} = 0;
        varargout{3} = errStr;
        return;
    end
    if ~isempty(information.DE_scheme)
        information.DE_scheme = squeeze(information.DE_scheme);
        
        if isfield(information,'qsi'),
            if information.qsi,
                information.DE_scheme = correctBvaluesForQSI(information.DE_scheme);
            end;
        end;

        
        [res, errStr]= dw_data_admin('set_diffDirs',information.DE_scheme); % write diff encod directions
        if ~isempty(errStr)
            if isempty(vHd)
                disp(errStr);
            else
                set(vHd, 'string', errStr);
                drawnow;
            end
            varargout{2} = 0;
            varargout{3} = errStr;
            return;
        end
    end
    if ~isempty(information.bfactor)
        [res, errStr]= dw_data_admin('set_bVal', information.bfactor); % write b-value, several b-values possible!
        if ~isempty(errStr)
            if isempty(vHd)
                disp(errStr);
            else
                set(vHd, 'string', errStr);
                drawnow;
            end
            varargout{2} = 0;
            varargout{3} = errStr;
            return;
        end
    end
    
    [res, errStr]= dw_data_admin('set_TR', information.TR); % write TR
    if ~isempty(errStr)
        if isempty(vHd)
            disp(errStr);
        else
            set(vHd, 'string', errStr);
            drawnow;
        end
        varargout{2} = 0;
        varargout{3} = errStr;
        return;
    end
    [res, errStr]= dw_data_admin('set_TE', information.TE); % write TE
    if ~isempty(errStr)
        if isempty(vHd)
            disp(errStr);
        else
            set(vHd, 'string', errStr);
            drawnow;
        end
        varargout{2} = 0;
        varargout{3} = errStr;
        return;
    end
    [res, errStr]= dw_data_admin('set_TI', information.TI); % write TI
    if ~isempty(errStr)
        if isempty(vHd)
            disp(errStr);
        else
            set(vHd, 'string', errStr);
            drawnow;
        end
        varargout{2} = 0;
        varargout{3} = errStr;
        return;
    end
    [res, errStr]= dw_data_admin('set_orient', information.PatientPosition); % ist das der richtige Eintrag?
    if ~isempty(errStr)
        if isempty(vHd)
            disp(errStr);
        else
            set(vHd, 'string', errStr);
            drawnow;
        end
        varargout{2} = 0;
        varargout{3} = errStr;
        return;
    end
    [res, errStr]= dw_data_admin('set_patient', information.patient); % write patient name
    if ~isempty(errStr)
        if isempty(vHd)
            disp(errStr);
        else
            set(vHd, 'string', errStr);
            drawnow;
        end
        varargout{2} = 0;
        varargout{3} = errStr;
        return;
    end
    [res, errStr]= dw_data_admin('set_vox', information.voxel); % write voxel size
    if ~isempty(errStr)
        if isempty(vHd)
            disp(errStr);
        else
            set(vHd, 'string', errStr);
            drawnow;
        end
        varargout{2} = 0;
        varargout{3} = errStr;
        return;
    end
    [res, errStr]= dw_data_admin('set_user', user); % write user struct with remaining information
    if ~isempty(errStr)
        if isempty(vHd)
            disp(errStr);
        else
            set(vHd, 'string', errStr);
            drawnow;
        end
        varargout{2} = 0;
        varargout{3} = errStr;
        return;
    end
    [res, errStr]= dw_data_admin('close');
    if ~isempty(errStr)
        if isempty(vHd)
            disp(errStr);
        else
            set(vHd, 'string', errStr);
            drawnow;
        end
        varargout{2} = 0;
        varargout{3} = errStr;
        return;
    else
        mes= 'info file saved';
        if isempty(vHd)
            disp(mes);
        else
            set(vHd, 'string', mes);
            drawnow;
        end
    end
else
    [res,errStr] = dw_data_admin('open',fName);
    if ~isempty(errStr)
        if isempty(vHd)
            disp(errStr);
        else
            set(vHd, 'string', errStr);
            drawnow;
        end
        varargout{2} = 0;
        varargout{3} = errStr;
        return;
    end
    user = dw_data_admin('get_user');
    information.DE_scheme = dw_data_admin('get_diffDirs');
    
    
    information.bfactor = dw_data_admin('get_bVal');
    [res, errStr]= dw_data_admin('close');
    if ~isempty(errStr)
        if isempty(vHd)
            disp(errStr);
        else
            set(vHd, 'string', errStr);
            drawnow;
        end
        varargout{2} = 0;
        varargout{3} = errStr;
        return;
    end
end

% function to read dicoms and put it into binary file
try
    if flag_averages == 0
        [res,errStr] = get_data_from_dicom(list,user.acquNo,user.SliceNo,...
            user.mosaic.status,user.mosaic.NumberOfImagesInRow,user.matrix(1),...
            user.matrix(2),fName,user.SliceLocation,information.DE_scheme,...
            information.edges,flag_averages,user.nob0s,information.DO_interpolation,[],vHd);
    else
        information.DO_interpolation = 0; % when averaging modus no interpolation is possible, this might cause problems for GE data, if somone wants to average
        [res,errStr] = get_data_from_dicom(list,user.acquNo,user.SliceNo,...
            user.mosaic.status,user.mosaic.NumberOfImagesInRow,user.matrix(1),...
            user.matrix(2),fName,user.SliceLocation,information.DE_scheme,...
            [],flag_averages,user.nob0s,information.DO_interpolation,count_averages,vHd);
    end
    if ~isempty(errStr)
        if isempty(vHd)
            disp(errStr);
        else
            set(vHd, 'string', errStr);
            drawnow;
        end
        varargout{2} = 0;
        varargout{3} = errStr;
        return;
    end
catch
    tmpStr= lasterr;
    tmpStr(tmpStr == char(10))= ' ';
    errStr = strcat('Error within function "read_DTI_Dicom" (get_data_from_dicom):',tmpStr);
    if isempty(vHd)
        disp(errStr);
    else
        set(vHd, 'string', errStr);
        drawnow;
    end
    varargout{2} = 0;
    varargout{3} = errStr;
    return;
end

mes= 'data file saved';
if isempty(vHd)
    disp(mes);
else
    set(vHd, 'string', mes);
    drawnow;
end
if isempty(information.DE_scheme)
    sts = 0;
else
    sts = 1;
end
if isempty(user.nob0s)
    sts(end+1) = 0;
else
    sts(end+1) = 1;
end
if isempty(information.bfactor)
    sts(end+1) = 0;
else
    sts(end+1) = 1;
end
if nargout == 1
    varargout{1} = fName; % name of output file
elseif nargout == 2
    varargout{1} = fName;
    varargout{2} = sts;
elseif nargout == 3
    varargout{1} = fName;
    varargout{2} = sts;
    varargout{3} = errStr;
elseif nargout > 3
    varargout{1} = fName;
    varargout{2} = sts;
    varargout{3} = errStr;
    for m = 4 : nargout
        varargout{m} = [];
    end
end

end
% end function read_DTI_Siemens
%_________________________________________________________________________
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%==========================================================================
% local functions

% function get_header_info
%--------------------------------------------------------------------------
    function [gather,user,fName] = get_header_info(list,vHd)

        user.nob0s = 0;
        user.SliceLocation = [];
        gather.bfactor = [];
        gather.DE_scheme = [];

        % read headers and extract necessary information for info file
        if isempty(vHd)
            wait = waitbar(0,'Please wait ... looking through dicoms.');
        end
        for m = 1: size(list,1)
            if isempty(vHd)
                waitbar(m/(size(list,1)))
            else
                msgStr= sprintf('Gather information from dicom file %d of %d.', m, size(list,1));
                set(vHd, 'String', msgStr);
                drawnow;
            end

            if exist(list(m,:).name,'file')~= 2
                error('Error in read_DTI_dicom: dicom file %s was not found.', list(m,:));
            end
            % read header and data of dicom file
            [DTI,header(m).data] = dicom_read_singlefile(list(m,:).name,1);
            gather.edges(:,:,m) = DTI.edges;

            % check if mosaic or not, seems to work for VA25, VB12 and VB13
            if m == 1
                if isfield(header(m).data,'Private_0029_1010')
                    for n=1:size(header(m).data.Private_0029_1010,2)
                        if strcmp(header(m).data.Private_0029_1010(n).name,'NumberOfImagesInMosaic')
                            InfoMosaic = header(m).data.Private_0029_1010(n);
                            break;
                        end
                    end
                    if isempty(InfoMosaic.item)
                        user.mosaic.status = 0;
                    else
                        user.mosaic.status = 1;
                    end
                else
                    user.mosaic.status = 0;
                end
                if strcmp(header(m).data.ProtocolName,'ek_dsi'),
                    gather.qsi = true;
                else
                    gather.false = true;
                end;
            end

            % find necessary information depending on mosaic or no mosaic
            if user.mosaic.status == 1
                user.SliceNo = str2double(InfoMosaic.item(1,1).val); % amount of slices from mosaic info header
                user.mosaic.NumberOfImagesInRow = fix(sqrt(user.SliceNo));
                while user.mosaic.NumberOfImagesInRow^2 < user.SliceNo
                    user.mosaic.NumberOfImagesInRow = user.mosaic.NumberOfImagesInRow +1;
                end
                Rows = header(m).data.Rows/user.mosaic.NumberOfImagesInRow;
                Columns = header(m).data.Columns/user.mosaic.NumberOfImagesInRow;
                user.acquNo = size(list,1);
                gather.DO_interpolation = 0;
                %user.InstanceNo(m) = [];
            else
                %determine sliceNo using SliceLocation
                if m == 1
                    user.SliceLocation = header.data.SliceLocation;
                    %count = 1;
                end
                % find all possible slice locations
                if m > 1 && ~any(user.SliceLocation == header(m).data.SliceLocation)
                    user.SliceLocation(end+1) = header(m).data.SliceLocation;
                end
%                 tmp_InstanceNo(m,:) = header(m).data.SOPInstanceUID;
%                 user.InstanceNo(m,:) = str2num(tmp_InstanceNo(m,end-7:end));
                user.acquNo = 0;
                user.SliceNo = 0;
                if isfield(header(m).data,'AcquisitionMatrix')
                    Rows = header(m).data.AcquisitionMatrix(1);
                    Columns = header(m).data.AcquisitionMatrix(4);
                    Rows_orig = header(1).data.Rows;
                    Columns_orig = header(1).data.Columns;
                    if Rows * 2 == Rows_orig && Columns * 2 == Columns_orig
                        gather.DO_interpolation = 1;
                        pixel = header(1).data.PixelSpacing.*2;
                    else
                        gather.DO_interpolation = 0;
                        Rows = header(1).data.Rows;
                        Columns = header(1).data.Columns;
                    end
                else
                    gather.DO_interpolation = 0;
                    Rows = header(m).data.Rows;
                    Columns = header(m).data.Columns;
                end
                user.mosaic.NumberOfImagesInRow = 1;
            end

            % extract Image Orientation (for two header versions)
            % dicom 3
            if isfield(header(m).data, 'ImageOrientationPatient')
                ImgPatOrient(:,m) = header(m).data.ImageOrientationPatient;
            end
            % ACR-NEMA_2.
            if isfield(header(m).data, 'ImageOrientation')
                ImgPatOrient(:,m) = header(m).data.ImageOrientation;
            end

            % extract diffusion sequence specific information, seems to work for
            % Siemens VA25 and VB12 and VB13!
            indDiffDir = []; infB_value = [];
            if isfield(header(m).data,'Private_0029_1010')
                for n = 1 : length(header(m).data.Private_0029_1010)
                    if strcmp(header(m).data.Private_0029_1010(n).name,'DiffusionGradientDirection')
                        indDiffDir = header(m).data.Private_0029_1010(n);
                    end
                    if strcmp(header(m).data.Private_0029_1010(n).name,'B_value' )
                        infB_value = header(m).data.Private_0029_1010(n);
                    end
                    if ~isempty(indDiffDir) && ~isempty(infB_value)
                        break;
                    end
                end
                if ~isempty(indDiffDir) && ~isempty(infB_value)
                    if infB_value.nitems > 0,
                        [gather.DE_scheme(1,m,:), bvalue] = extract_DWI_parameters(header(m).data, DTI.edges);
                        
                        gather.DE_scheme(1,m,:) = gather.DE_scheme(1,m,:)*sqrt(bvalue); % to save bval information (mrc)
                        
                        if bvalue == 0
                            user.nob0s = user.nob0s + 1;
                        elseif indDiffDir.nitems > 0
                            if ~isempty(gather.bfactor) && gather.bfactor(end) ~= bvalue
                                gather.bfactor(end+1) = bvalue;
                            else
                                gather.bfactor = bvalue;
                            end
                        else
                            gather.DE_scheme = [];
                            gather.bfactor = [];
                            user.nob0s = [];
                        end
                    end;
                else
                    gather.DE_scheme = [];
                    gather.bfactor = [];
                    user.nob0s = [];
                end
            else
                gather.DE_scheme = [];
                gather.bfactor = [];
                user.nob0s = [];
            end

            % calculate Acquistion Time
            if m == 1                                                                % needs to be tested if really always first image!!
                hour1 = str2double(header(m).data.AcquisitionTime(1,1:2))*60*60;
            end
            hours = (str2double(header(m).data.AcquisitionTime(1,1:2))*60*60-hour1); % in seconds % if in milliseconds: * 1000
            minutes = str2double(header(m).data.AcquisitionTime(1,3:4))*60;          % in seconds % if in milliseconds: * 1000
            seconds = str2double(header(m).data.AcquisitionTime(1,5:6));             % in seconds % if in milliseconds: * 1000
            micro = str2double(header(m).data.AcquisitionTime(1,8:end))/1000000;     % in seconds % if in milliseconds: /1000
            user.AcquisitionTime(m) =  hours + minutes + seconds + micro;            % in seconds, each image has its own time!
        end % for
        if isempty(vHd)
            close(wait)
        end

        if user.mosaic.status == 0 % no mosaic
            user.SliceNo = size(user.SliceLocation,2); % finally the number of slices
            user.acquNo = size(list,1)/user.SliceNo; % resulting in number of acquisitons as well
            user.SliceLocation = sort(user.SliceLocation); % sort Slice Location in order to have Slice ID
            %user.InstanceNo = sort(user.InstanceNo); % InstanceNo only for
            %testing if sorting is correct
        end

        % for no mosaic case and existing DE information save only DE_scheme for one whole slice stack
        if user.mosaic.status == 0 && ~isempty(gather.DE_scheme)
            ind = [1 : user.SliceNo : size(list,1)];
            gather.DE_scheme = gather.DE_scheme(1,ind,:);
            user.nob0s = user.nob0s/user.SliceNo;
        end

        % Patients Name, TI, TR, TE etc for info file!
        if ischar(header(1).data.PatientsName)
            gather.patient = header(1).data.PatientsName;
        else
            gather.patient = [header(1).data.PatientsName.GivenName,' ', header(1).data.PatientsName.FamilyName];
        end
        if isfield(header(1).data, 'InversionTime')
            gather.TI = header(1).data.InversionTime;
        else
            gather.TI = 0;
        end
        
        if isfield(header(1).data, 'RepetitionTime')
            gather.TR = header(1).data.RepetitionTime;
        else
            gather.TR = 0;
        end
        if isfield(header(1).data, 'EchoTime')
            gather.TE = header(1).data.EchoTime;
        else
            gather.TE = 0;
        end
        if isfield(header(1).data, 'PatientPosition')
            gather.PatientPosition = header(1).data.PatientPosition;
        else
            gather.PatientPosition = 'none';
        end     
        if user.mosaic.status == 0 && gather.DO_interpolation == 1
            gather.voxel = [pixel', header(1).data.SliceThickness];
        else
            gather.voxel = [header(1).data.PixelSpacing', header(1).data.SliceThickness];
        end
        
        % information from header which should be the same for all Versions (VA25,
        % VB12 and VB13)
        if isfield(header(1).data, 'FlipAngle')
            user.FlipAngle = header(1).data.FlipAngle;
        else
            user.FlipAngle = 0;
        end
        
        if user.mosaic.status == 1
            user.ImagePatientOrientation = reshape(ImgPatOrient,3,2,user.acquNo);
            gather.edges = DTI.edges;
        else
            user.ImagePatientOrientation = reshape(ImgPatOrient,3,2,size(list,1));
        end
        if isfield(header(1).data, 'MagneticFieldStrength')
            user.MagneticFieldStrength = header(1).data.MagneticFieldStrength;
        else
            user.MagneticFieldStrength = 0;
        end
        if isfield(header(1).data, 'Manufacturer')
            user.Manufacturer = header(1).data.Manufacturer;
        else
            user.Manufacturer = 0;
        end
        if isfield(header(1).data, 'StudyDate')
            user.ScanDate = header(1).data.StudyDate;
        else
            user.ScanDate = 0;
        end
        if isfield(header(1).data, 'ProtocolName')
            user.SequenceName = header(1).data.ProtocolName;
        else
            user.SequenceName = 0;
        end
        if isfield(header(1).data, 'SoftwareVersions')
            user.SoftwareVersion = header(1).data.SoftwareVersions;
        else
            user.SoftwareVersion = 0;
        end

        user.matrix = [Rows,Columns];

        % create name of info file and binary data file
        index = strfind(list(1,1).name,'\');
        patient = strcat(gather.patient,'_',num2str(user.ScanDate));
        ind_points = strfind(patient,'.');
        if ~isempty(ind_points)
            patient(ind_points) = [];
        end
        ind_commas = strfind(patient,',');
        if ~isempty(ind_commas)
            patient(ind_commas) = [];
        end
        patient(find(isspace(patient)==1)) = [];
        if isempty(index)
            index = strfind(list(1,1).name,'/');
            if isempty(index)
                fName = 'test';
            else
                fName = strcat(list(1,1).name(1:index(end)),patient);
            end
        else
            fName = strcat(list(1,1).name(1:index(end)),patient);
        end
        
    end
% end function get_header_info
%--------------------------------------------------------------------------


% funtion get_data_from_dicom
%--------------------------------------------------------------------------
    function [res,errStr] = get_data_from_dicom(list,acquNo,sliceNo,mosaic,NumberOfImagesInRow,Rows,Columns,fName,SliceLocation,DE_scheme,edges_temp,flag_averages,nob0s,DO_interpolation,average_count,vHd)
    
    [res, errStr]= dw_data_admin('open',fName);
    if ~isempty(errStr)
        if isempty(vHd)
            disp(errStr);
        else
            set(vHd, 'string', errStr);
            drawnow;
        end
        res = -1;
        return
    end

    if flag_averages == 0
        edges = zeros(size(edges_temp));
    end
    if isempty(vHd)
        if average_count > 1
            waitstring = ['Please wait ... write average ',num2str(average_count),' into binary file.'];
            wait = waitbar(0,waitstring);
        else
            wait = waitbar(0,'Please wait ... writing dicom files into binary file.');
        end
    end
    % save DTI data in file
    if mosaic == 0              % no mosaic images;
        count_slices = sliceNo;
        count_b0 = 0;
        for m = 1 : size(list,1)
            [DTI,header,errStr] = dicom_read_singlefile(list(m,:).name,1);
            if ~isempty(errStr)
                if isempty(vHd)
                    disp(errStr);
                else
                    set(vHd, 'string', errStr);
                    drawnow;
                end
                res = -1;
                return
            end
           % pixel interpolation mostly for GE data
            if DO_interpolation == 1
                tmp = interpft(DTI.dataAy, Rows , 1);
                DTI.dataAy = interpft(tmp ,Columns , 2);
            end
            % put data to correct slice position
            slice_id = find(SliceLocation == header.SliceLocation);
            % find number of acquisition and put data to correct place
            if ~isempty(DE_scheme)
                if isfield(header, 'ImageOrientationPatient')
                    ImgPatOrient = header.ImageOrientationPatient;
                end
                % ACR-NEMA_2.
                if isfield(header, 'ImageOrientation')
                    ImgPatOrient = header.ImageOrientation;
                end
                [DEdir,b_value,errStr] = extract_DWI_parameters(header, DTI.edges);
                if ~isempty(errStr)
                    if isempty(vHd)
                        disp(errStr);
                    else
                        set(vHd, 'string', errStr);
                        drawnow;
                    end
                    res = -1;
                    return;
                end
                for p = 1 : size(DE_scheme,1)
                    % if more than one b0 image this has to be taken in account
                    % when looking for correct acquNo, not tested yet SuS 08.02.2008
                    if b_value == 0 && nob0s > 1
                        count_slices = count_slices - 1;
                    end
                    if count_slices == 0 && b_value == 0
                        count_b0 = count_b0 + 1;
                        if DE_scheme(p,:) == DEdir
                            acquNo_id(end+1) = p;
                        end
                        acquNo_id = acquNo_id(count_b0);
                        if count_b0 == nob0s
                            break;
                        end
                    else
                        if DE_scheme(p,:) == DEdir
                            acquNo_id = p;
                            break;
                        end
                    end
                end
                % add data if flag_averages equals 1
                if flag_averages == 1
                    [res, errStr]= dw_data_admin('addData', slice_id, acquNo_id ,DTI.dataAy);
                else
                    [res, errStr]= dw_data_admin('putData', slice_id, acquNo_id ,DTI.dataAy);
                end
                if ~isempty(errStr)
                    if isempty(vHd)
                        disp(errStr);
                    else
                        set(vHd, 'string', errStr);
                        drawnow;
                    end
                    res = -1;
                    return;
                end
                if flag_averages == 0
                    edges(:,:,slice_id+(sliceNo*(acquNo_id-1))) = edges_temp(:,:,m);
                end
            else % others than Siemens
                % SuS Feb 2008: no averaging from seperate file folders possible yet, is this necessary at all?
                % find out order of data and sort into binary file
                if m == 1
                    acquNo_id = 1;
                    [res, errStr]= dw_data_admin('putData', slice_id, acquNo_id,DTI.dataAy);
                elseif m > 1 && slice_id_old == slice_id % all diffDirs were colected for one slice first
                    acquNo_id = acquNo_id + 1;
                    [res, errStr]= dw_data_admin('putData', slice_id, acquNo_id,DTI.dataAy);
                    if acquNo_id == acquNo
                        acquNo_id = 0;
                    end
                elseif m > 1 && slice_id_old ~= slice_id % one diffDir was collected for all slices first
                    [res, errStr]= dw_data_admin('putData', slice_id, acquNo_id,DTI.dataAy);
                    if slice_id == sliceNo
                        acquNo_id = acquNo_id + 1;
                    end
                end
                slice_id_old = slice_id;
                edges = edges_temp; % SuS Feb 2008: think more about edges, here just the order of files is taken...
            end
            if ~isempty(errStr)
                if isempty(vHd)
                    disp(errStr);
                else
                    set(vHd, 'string', errStr);
                    drawnow;
                end
                res = -1;
                return
            end
            if flag_averages == 1
                msgStr= sprintf('Added dicom image %d/%d to binary file, average file %d.', m, size(list,1),average_count);
            else
                msgStr= sprintf('Saved dicom image %d/%d in binary file.', m, size(list,1));
            end
            if isempty(vHd)
                waitbar(m/(size(list,1)))
            else
                set(vHd, 'String', msgStr);
                drawnow;
            end
        end % for
        if flag_averages == 0
            edges = get_volume_geo(edges);
        end
    else                        % mosaic images
        for m = 1: size(list,1)
            if flag_averages == 1
                outStr= sprintf('added DTI data (%d/%d) to binary file, average file %d.', m, acquNo, average_count);
            else
                outStr= sprintf('save DTI data in binary file ... (%d/%d)', m, acquNo);
            end
            if isempty(vHd)
                waitbar(m/(size(list,1)))
            else
                set(vHd, 'string', outStr);
                drawnow;
            end
            [DTI_volume, header, errStr] = dicom_read_singlefile(list(m,:).name, 1);
            if ~isempty(errStr)
                if isempty(vHd)
                    disp(errStr);
                else
                    set(vHd, 'string',errStr);
                    drawnow;
                end
                res = -1;
                return;
            end
            % find correct position of slices in mosaic
            for Slice = 1 : sliceNo
                RowNumber= fix(Slice/NumberOfImagesInRow - 0.00001) + 1;
                ColNumber = Slice - (RowNumber - 1)* NumberOfImagesInRow;
                RowStart = 1 + Rows * (RowNumber-1);
                ColStart = 1 + Columns * (ColNumber-1);
                volume(:,:,Slice) = DTI_volume.dataAy(RowStart : RowStart + Rows -1 ,ColStart : ColStart + Columns -1, : );
            end
            % write DTI volume in binary file
            if flag_averages == 1
                [res, errStr]= dw_data_admin('addData', 1:sliceNo, m ,volume);
            else
                [res, errStr]= dw_data_admin('putData', 1:sliceNo, m ,volume);
            end
            if ~isempty(errStr)
                if isempty(vHd)
                    disp(errStr);
                else
                    set(vHd, 'string', errStr);
                    drawnow;
                end
                res = -1;
                return;
            end
        end
        edges = DTI_volume.edges;
    end
    if isempty(vHd)
        close(wait)
    end
    
    % flag_averages is also 0, if the first of the data sets to average was
    % only processed yet
    if flag_averages == 0
        [res, errStr]= dw_data_admin('set_edges', edges); % write edges derived from dicom_read_singlefile
        if ~isempty(errStr)
            if isempty(vHd)
                disp(errStr);
            else
                set(vHd, 'string', errStr);
                drawnow;
            end
            res = -1;
            disp('Error in read_dti_dicom: Could not write edges, either not availible from header or wrong format');
            return;
        end
    end
    [res, errStr]= dw_data_admin('close');
    if ~isempty(errStr)
        if isempty(vHd)
            disp(errStr);
        else
            set(vHd, 'string', errStr);
            drawnow;
        end
        res = -1;
        return;
    else
        res = 1;
        errStr = [];
    end
    end
    % end function get_data_from_dicom
    %--------------------------------------------------------------------------


% function extract_DWI_parameters
%--------------------------------------------------------------------------
    function [DEdir,b_value,errStr] = extract_DWI_parameters(infoMx, edges)

        % extract DiffDirection and b-value
        infoDiff = [];
        infoBvalue = [];
        for m = 1 : length(infoMx.Private_0029_1010)
            if strcmp(infoMx.Private_0029_1010(m).name,'DiffusionGradientDirection')
                infoDiff = infoMx.Private_0029_1010(m);
            end
            if strcmp(infoMx.Private_0029_1010(m).name,'B_value' )
                infoBvalue = infoMx.Private_0029_1010(m);
            end
            if ~isempty(infoDiff) && ~isempty(infoBvalue)
                break;
            end
        end
        

        % b value
        if isempty(infoBvalue.item)
            b_value = [];          % really necessary????
            errStr = 'Error in extract_DWI_parameters: no entry in Private_002910 with field "B_value"';
            if isempty(infoDiff.item)
               %errStr = 'Error in extract_DWI_parameters: no diffusion information in header!'; 
               DEdir = [];
               errStr = 'Error in extract_DWI_parameters: no entry in Private_002910 with field "B_value" and "DiffusionGradientDirection"';
            end
            return;
        else
            b_value = str2double(infoBvalue.item(1).val);
        end

        %DE direction
        if b_value == 0
            DEdir = [0 0 0];
        elseif ~isempty(infoDiff.item)
            DEdir = str2num([infoDiff.item(1).val,infoDiff.item(2).val,infoDiff.item(3).val]);
            % rotation of Diffusion endcoding vector in image coordinates
            if isempty(edges)
                errStr = 'Error in extract_DWI_parameters: parameter "edges" missing.';
                return;
            end
            [p_org,errStr]= hm_analyse(edges);
            if ~isempty(errStr)
                b_value = [];
                DEdir = [];
                disp(errStr);
                return
            end
            p_rot= [0 0 0 p_org(4:6) sign(p_org(7:9)) 0 0 0];
            rotHM= hm_create(p_rot);
            rotM= rotHM(1:3, 1:3)';
            DEdir = rotM*DEdir';
            DEdir = DEdir';
        end
        errStr = [];
    end

% end function extract_DWI_parameters
%--------------------------------------------------------------------------










function newDE_scheme = correctBvaluesForQSI(DE_scheme)

    h = load('allhex.mat');

    schemelen = cellfun(@(x) length(x.bfactor)-1,h.allhex);

    b = squeeze(sum(DE_scheme.^2,2));
    bidx = find(b>0);

    hexschidx = find(schemelen==length(bidx));
    hexb = h.allhex{hexschidx}.bfactor;

    fprintf('QSI scheme detected!\n number of q-points :%i\n ',length(bidx));

    newDE_scheme = DE_scheme;
    for k = 1:length(bidx)
         dir = DE_scheme(bidx(k),:);
         dir = dir*sqrt(hexb(k+1));
         newDE_scheme(bidx(k),:) = dir;
    end;

end




% end local functions
%==========================================================================


