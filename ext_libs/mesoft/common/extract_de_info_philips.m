% [gather user] = extract_de_info_philips
%

% Susanne Schnell 17.03.2009

function [gather user] = extract_de_info_philips

path = uigetdir(pwd,'Select the directory where the dicom files are located.'); % Matlab 2007b
%path = 'E:\MR\DTI\DataSets\ThimoGrotz_15012008\DTI_slices69_vox2x2x2_matr104x104_TR11000_TE94_DEdirs61_8b0s'; % Matlab R14
if isequal(path , 0)
    frprintf('User pressed cancel');
    return;
else
    fprintf('User selected dcm files from %s.', path);
end

if ispc
    files = dir(strcat(path,'\*.dcm'));
else
    files = dir(strcat(path,'/*.dcm'));
end

user.nob0s = 0;
user.SliceLocation = [];
gather.bfactor = [];

% read headers and extract necessary information for info file
wait = waitbar(0,'Please wait ... looking through dicoms.');

for m = 1: size(files,1)
    waitbar(m/(size(files,1)))

    if exist(files(m,:).name,'file')~= 2
        error('Error in read_dti_dicom: dicom file %s was not found.', files(m,:));
    end
    % read header and data of dicom file
    [DTI,header(m).data] = dicom_read_singlefile(files(m,:).name,1);
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
        user.acquNo = size(files,1);
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
    if isfield(header(m).data,'Private_2005_10B0')
        gather.DE_scheme(1,m,1) = header(m).data.Private_2005_10B0;
        gather.DE_scheme(1,m,2) = header(m).data.Private_2005_10B1;
        gather.DE_scheme(1,m,3) = header(m).data.Private_2005_10B2;
        gather.bfactor(m) = header(m).data.Private_2001_1003;

        % rotate DEdirs
        [p_org,errStr]= hm_analyse(gather.edges(:,:,m));
        p_rot= [0 0 0 p_org(4:6) sign(p_org(7:9)) 0 0 0];
        rotHM= hm_create(p_rot);
        rotM= rotHM(1:3, 1:3)';
        DEdir = rotM * squeeze(gather.DE_scheme(1,m,:));
        gather.DE_scheme(1,m,:) = DEdir';
    end

    if gather.bfactor(m) == 0;
        user.nob0s = user.nob0s + 1;
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

close(wait)

if user.mosaic.status == 0 % no mosaic
    user.SliceNo = size(user.SliceLocation,2); % finally the number of slices
    user.acquNo = size(files,1)/user.SliceNo; % resulting in number of acquisitons as well
    user.SliceLocation = sort(user.SliceLocation); % sort Slice Location in order to have Slice ID
    %user.InstanceNo = sort(user.InstanceNo); % InstanceNo only for
    %testing if sorting is correct
end

% for no mosaic case and existing DE information save only DE_scheme for one whole slice stack
if user.mosaic.status == 0 && ~isempty(gather.DE_scheme)
    ind = [1 : user.SliceNo : size(files,1)];
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
    user.ImagePatientOrientation = reshape(ImgPatOrient,3,2,size(files,1));
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
DE_scheme = squeeze(gather.DE_scheme);
save(fullfile(path,'DE_scheme.mat'),'DE_scheme');

        
        
  