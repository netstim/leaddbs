function [] = uw_overlay_convert2dicom(dicom_file, merged_file, newSeriesNumber, newSeriesDescription, outputDirectory, mergedImageVolume, outputImagePosition, outputOrientation) %#ok<INUSD>

%   FUNCTION uw_overlay_convert2dicom
%   Version 1.1
%
%   A function to take a merged NIfTI or ANALYZE dataset (anatomical + overlay
%   info) and convert it back to DICOM in axial LPS 'radiological orientation'.
%   The user specifies an original DICOM file, which is used to obtain patient,
%   study, and position information, and convert2dicom will create a new image
%   series with the merged dataset.
%
%   WARNING: I strongly suggest you read the ENTIRE documentation
%   (especially the notes on outputImagePosition) before you use it!
%
%   DEPENDANCIES: MATLAB NIfTI Tools (Included): load_nii, load_nii_hdr, xform_nii
%                 Does not require SPM.
%
%   SYNTAX:
%   [] = uw_overlay_convert2dicom(dicom_file, merged_file, newSeriesNumber,
%                  newSeriesDescription, outputDirectory, mergedImageVolume,
%                  outputImagePosition, outputOrientation)
%
%
%   *First three arguments are required, others optional.
%                  
%   dicom_file - String: Original DICOM file from the study you are using.
%                The patient name, study description, date, etc. will be
%                copied from this file.  Always specify the first file in a
%                series of images.  If possible, specify the same
%                anatomical series used to create the merged image.
%
%   merged_file - String: A NIfTI or ANALYZE image
%
%   newSeriesNumber - Number: Specify a NEW series number for the merged
%                     images.  Its very important that this be a UNIQUE
%                     number, or it will 'collide' with a previous series
%                     in your dataset and corrupt the AW database. And
%                     then your IT person won't like you anymore :-(
%
%   newSeriesDescription - String, optional.  Give the overlay description
%                          (such as 'fingertap left hand', etc).  Defaults
%                          to 'Merged Overlay' if not specified. 
%
%   outputDirectory - String, optional.  Specify the directory you would like
%                      to output the DICOM files.  Defaults to
%                      newSeriesNumber if not specified.  The DICOM files
%                      themselves will always be saved as I0001.dcm in
%                      ascending numerical order.
%
%   mergedImageVolume - Number, optional.  If merged_file is actually a
%                       series of image volumes (such as a time series),
%                       specify which volume you would like to use.
%                       Default is 1.
%               
%   outputImagePosition - Last but not least!  An optional argument to
%                         specify where convert2dicom should get its
%                         information regarding voxel size and the location
%                         of the origin.  Most of this is to deal with
%                         older or non-compliant programs, 0 should work fine
%                         with current versions of AFNI or SPM.
%
%           0 [default] - use the voxel size and origin location from the
%                         NIfTI/ANALYZE merged image IF it is marked with
%                         sform_code = 1 (NIFTI_XFORM_SCANNER_ANAT), which
%                         indicates scanner-based coordinate system.  If
%                         the data is not marked this way, the user
%                         is prompted if (s)he would like to use the
%                         NIfTI/ANALYZE header or the DICOM header for the
%                         positon information (or abort and think it over).
%
%                         *This setting isn't good for automation, since it
%                         may require user input.
%
%           1 - force NIfTI/ANALYZE position: Will always use the origin
%               and voxel size specified in the merged_file header.
%               Use this setting if you're sure that the NIfTI contains the
%               correct origin location and voxel size. This would
%               typically be true as long as sform_code = 1, and as long as
%               image was not realigned or modified after being imported
%               from the original anatomical DICOM.
%
%           2 - force DICOM position: Will always use the origin and voxel
%               size specified in the DICOM header.  Use this option to
%               replace the image data in the DICOM with NIfTI data,
%               without touching any of the position information.  It will
%               write the data in LPS orientation, regardless of the
%               orientation specified in the DICOM header
%
%          > 3 - NOT USED-> Future versions will allow more flexibility to
%                           decide how the position information will be handled,
%                           and will allow the user to override values by hand.
%
%   outputOrientation - NOT USED-> In future versions of this function, you
%                       will be able to specify the orientation you want to
%                       output the image, or have convert2dicom match the output
%                       image orientation with the original orientation of the
%                       DICOM.  Right now, all images are converted to axial LPS
%                       'radiological orientation', which is the default
%                       coordinate system used in DICOM:
%                       LPS = X: increases right     --> left
%                             Y: increases anterior --> posterior
%                             Z: increases inferior  --> superior
%
%        *For more information on orientations and coordinate systems see:
%         http://lcni.uoregon.edu/~jolinda/MRIConvert/fileformats.htm
%
%
%   Copyright 2007 Samuel A. Hurley - samuel.hurley[at]gmail.com
%   Universtiy of Wisconsin - Madison | Applied Neuro fMRI Lab
%   Published under GNU General Public License Version 3, see <http://www.gnu.org/licenses/>


%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License version 2 or 3 ONLY
%     as published by the Free Software Foundation
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% << Actual code below >>

% I. User input error handling.

% Check for the correct number of arguments
if nargin < 3       % First three arguments are required
    error('Arguments: You must specify at least three arguments.');
end

if nargin > 8       % Cannot have more than 8 args
    error('Arguments: You cannot specify more than eight arguments.  Check syntax.');
end

%Check that the first two arguments are strings & that the files exist:
if ~ ischar(dicom_file)
    error('Arguments: dicom_file must be a filename (string)');
end

if exist(dicom_file) ~= 2     %#ok<EXIST>  % exist returns 2 for a valid filename
    error(['Arguments:' dicom_file ' - file not found.']);
end

if ~ ischar(merged_file)
    error('Arguments: merged_file must be a filename (string)');
end

if exist(merged_file) ~= 2  %#ok<EXIST>
    error(['Arguments:' merged_file ' - file not found.']);
end

%Check that newSeriesNumber is indeed a number
if ~ isnumeric(newSeriesNumber)
    error('Arguments: newSeriesNumber must be a number.');
end

%Check that its a positive integer
if newSeriesNumber < 1 || ~ check_integer(newSeriesNumber)
    error('Arguments: newSeriesNumber must be a positive integer.');
end

% Check if optional arguments are specified, otherwise assign default values.

switch nargin
    case 3          % No optional args, use defaults
        newSeriesDescription = 'Merged Overlay';
        outputDirectory = num2str(newSeriesNumber);
        mergedImageVolume = 1;
        outputImagePosition = 0;
        outputOrientation = 'LPS';  %#ok<NASGU> % Not currently used, reserved for future versions
    case 4
        if ~ ischar(newSeriesDescription)
            error('Arguments: newSeriesDescription must be a string');
        end
        outputDirectory = num2str(newSeriesNumber);
        mergedImageVolume = 1;
        outputImagePosition = 0;
        outputOrientation = 'LPS';  %#ok<NASGU> % Not currently used, reserved for future versions
    case 5
        if ~ ischar(newSeriesDescription)
            error('Arguments: newSeriesDescription must be a string');
        end
        if ~ ischar(outputDirectory)
            error('Arguments: outputDirectory must be a string');
        end
        mergedImageVolume = 1;
        outputImagePosition = 0;
        outputOrientation = 'LPS';  %#ok<NASGU> % Not currently used, reserved for future versions
    case 6
        if ~ ischar(newSeriesDescription)
            error('Arguments: newSeriesDescription must be a string');
        end
        if ~ ischar(outputDirectory)
            error('Arguments: outputDirectory must be a string');
        end
        % mergedImageVolume must be a positive integer:
        if mergedImageVolume < 1 || ~ check_integer(mergedImageVolume)
            error('Arguments: mergedImageVolume must be a positive integer.');
        end
        outputImagePosition = 0; 
        outputOrientation = 'LPS';   %#ok<NASGU> % Not currently used, reserved for future versions
    case 7
        if ~ ischar(newSeriesDescription)
            error('Arguments: newSeriesDescription must be a string');
        end
        if ~ ischar(outputDirectory)
            error('Arguments: outputDirectory must be a string');
        end
        % mergedImageVolume must be a positive integer:
        if mergedImageVolume < 1 || ~ check_integer(mergedImageVolume)
            error('Arguments: mergedImageVolume must be a positive integer.');
        end
        % outputImagePosition must be 0, 1, or 2
        if ~ max(outputImagePosition == [0 1 2])
            error('Arguments: outputImagePosition must be a numerical value of 0, 1, or 2');
        end
        outputOrientation = 'LPS';  %#ok<NASGU> % Not currently used, reserved for future versions
    case 8
        if ~ ischar(newSeriesDescription)
            error('Arguments: newSeriesDescription must be a string');
        end
        if ~ ischar(outputDirectory)
            error('Arguments: outputDirectory must be a string');
        end
        % mergedImageVolume must be a positive integer:
        if mergedImageVolume < 1 || ~ check_integer(mergedImageVolume)
            error('Arguments: mergedImageVolume must be a positive integer.');
        end
        % outputImagePosition must be 0, 1, or 2
        if max(outputImagePosition == [0 1 2])
            error('Arguments: outputImagePosition must be a numerical value of 0, 1, or 2');
        end
        % Warn the user that outputOrientation dosen't actually do anything, yet.
        disp('Warning: Arguemnts: outputOrientation is not yet supported.  Setting this option will');
        disp('         not actually do anything.  Output orientation will always be LAS (radiological).');
end

% II. Declare variables, load stuff

% There are a lot of variables associated with position, voxel size, and
% orientation from the DICOM and NIfTI input files.  So, we will use the
% datastructure imageParams to group all of the values read from each file

imageParamsDICOM = struct('Orientation', [], ...
                          'ImageDim', [], ...       % Size of image array
                          'VoxelDim', [], ...       % Size of voxels
                          'Origin', []   ...       % The xyz of lower-right voxel
                          );
                     
imageParamsNIfTI = struct('Orientation', [1 0 0 0 1 0], ...  % Will always be converted to LAS
                          'ImageDim', [], ...
                          'VoxelDim', [], ...
                          'Origin', []   ...
                           );

% Load header data from DICOM and NIfTI
imageDICOM = dicominfo(dicom_file);
imageNIfTI = ea_load_nii(merged_file); % load_nii, from the MATLAB NIfTI tools software
                                    % will load any ANALYZE or NIfTI image
                                    % and convert it to RAS coordinates,
                                    % preserving the origin and voxel size
                                    % information.

% Calls function to fill out the imageParams structures, displays warnings
% if things are not loaded correctly:
imageParamsDICOM = getDicomParams(imageDICOM, imageParamsDICOM);
imageParamsNIfTI = getNiftiParams(imageNIfTI, imageParamsNIfTI);


% III: Display some information from the DICOM and NIfTI to the user, pass on warnings
show_userinfo(imageDICOM, imageParamsDICOM, imageNIfTI, imageParamsNIfTI);



% IV: Check the value of outputImagePosition to determine wheter to use the
%     position and orientation information from the NIfTI header, the DICOM
%     header, or prompt the user

% Some temporary switches to set if the user wishes to override the NIfTI
% header's position and orientation information with those specified in the
% DICOM header:
useDICOMPosition = false(1);
useDICOMOrientation = false(1);

% Check the user's outputImagePosition preference:
switch outputImagePosition
    case 0 % Use the NIfTI header, as long as its marked as scanner coordinates
        % If sform is not zero, xform_nii will reorient using the sform matrix:
        if imageNIfTI.original.hdr.hist.sform_code > 0
            if imageNIfTI.original.hdr.hist.qform_code ~= 1 % Check to see if the code is SCANNER_ANATOMICAL
                form_code = true; % True means it ISNT defined correctly
            else
                form_code = false; % FALSE means it IS defined correctly.  Sorry, its confusing :o)
            end
        % Otherwise xform_nii will reorient using the qform matrix:
        elseif imageNIfTI.original.hdr.hist.sform_code > 0
            if imageNIfTI.original.hdr.hist.sform_code ~= 1
                form_code = true;
            else
                form_code = false;
            end
        else
            form_code = false;
        end
        
        % Prompt the user if there is a problem, otherwise set useDICOMPosition and Orientation to false
        if imageNIfTI.filetype == 0 || form_code
            disp('WARNING: There are one or more problems with your NIfTI header, ');
            disp('         please review the messages above for further information.');
            
            % Prompt the user which orientation to use
            promptOrientation = 'X';
            while  ~isequal(lower(promptOrientation), 'n') && ~isequal(lower(promptOrientation), 'd') && ~isequal(lower(promptOrientation), 'a')
                disp('Which orientation should I use?  N = NIfTI Header, D = DICOM Header, A = Abort');
                promptOrientation = input('[n/d/a]: ', 's');
            end

            if isequal(lower(promptOrientation), 'n')
                useDICOMOrientation = false(1);
            elseif isequal(lower(promptOrientation), 'd')
                useDICOMOrientation = true(1);
            else
                error('Error: Aborted by user');
            end
            
            promptPosition = 0;
            while ~isequal(lower(promptPosition), 'n') && ~isequal(lower(promptPosition), 'd') && ~isequal(lower(promptPosition), 'a')
                disp('Which position (origin+voxel size) should I use?  N = NIfTI Header, D = DICOM Header, A = Abort');
                promptPosition = input('[n/d/a]: ', 's');
            end

            if isequal(lower(promptPosition), 'n')
                useDICOMPosition = false(1);
            elseif isequal(lower(promptPosition), 'd')
                useDICOMPosition = true(1);
            else
                error('Error: Aborted by user');
            end
            
        else  % Everything is fine with the NIfTI header, use the NIfTI values:
            useDICOMPosition = false(1);
            useDICOMOrientation = false(1);
        end
        
    case 1 % Force use of the NIfTI header, warn the user if there might be a problem
        useDICOMPosition = false(1);
        useDICOMOrientation = false(1);
    case 2 % Force use of the DICOM header... this will just leave the DICOM header as it is 
           % and merely replace the image data with the NIfTI data instead.
        useDICOMPosition = true(1);
        useDICOMOrientation = true(1);
        % Warn the user if the DICOM isn't in LPS orientation:
        if ~ isequal(imageParamsDICOM.Orientation,[1 0 0 0 1 0])
            disp('WARNING: DICOM orientation is NOT LPS, but the new image data will be written in LPS orientation');
            disp('         This may mean that the DICOM you specified does not contain the original anatomical reference image');
            disp('         used to create your merged overlay.');    
        end
end


% V: Write a new DICOM header based on the user's preferences, and write the
%    image data into a series of new DICOM files.

% At this point we no longer need to pull information from the imageDICOM
% header, so we can start overwriting our new values into the variable.

% Not working yet... MATLAB's DICOMWRITE function refuses to write these
% two fields.  May need to modify the builtin code to support this option
% ########## Change DICOM implementation here ############################
% ########################################################################
% Change these two settings to match the implementation class used by your
% scanner system or PACS network:
  imageDICOM.ImplementationClassUID    = '1.2.840.113619.6.144';
  imageDICOM.ImplementationVersionName = 'AW4_2_06_02_EXT';
% ########################################################################

% Update the DICOM header with the new series description and UID
imageDICOM.SeriesNumber      = newSeriesNumber;
imageDICOM.SeriesDescription = newSeriesDescription;
imageDICOM.SeriesInstanceUID = dicomuid;  % Generate a unique identifier for the image series

% We always replace the DICOM rows/columns/slices with the NIfTI
% parameters, since these will match the dimensions of the new image
imageDICOM.Rows    = imageParamsNIfTI.ImageDim(1);
imageDICOM.Columns = imageParamsNIfTI.ImageDim(2);
imageDICOM.ImagesInAcquisition = imageParamsNIfTI.ImageDim(3);

% Write in the new position & orientation, if the user desires:
if useDICOMPosition == false
	imageDICOM.PixelSpacing(1)    = imageParamsNIfTI.VoxelDim(1);
	imageDICOM.PixelSpacing(2)    = imageParamsNIfTI.VoxelDim(2);
    imageDICOM.SliceThickness       = imageParamsNIfTI.VoxelDim(3);
    imageDICOM.SpacingBetweenSlices = imageParamsNIfTI.VoxelDim(3);
    
    imageDICOM.ImagePositionPatient = imageParamsNIfTI.Origin;
    imageDICOM.SliceLocation        = imageParamsNIfTI.Origin(3);
end

if useDICOMOrientation == false
    imageDICOM.ImageOrientationPatient = imageParamsNIfTI.Orientation';
end

% Create output directory and switch to it:
mkdir(outputDirectory);

% Grab the specified mergedImageVolume into the variable outputImageMatrix:
if length(size(imageNIfTI.img)) == 4
    outputImageMatrix = imageNIfTI.img(:,:,:,mergedImageVolume);
else
    outputImageMatrix = imageNIfTI.img;     % The matrix representing the output image
end

numberOfSlices = size(outputImageMatrix);   % Figure out the number of slices in the volume
numberOfSlices = numberOfSlices(3);

if numberOfSlices > 9999    % The DICOM filename cannot support more than 9999 images
    error('Error: Too many slices in the image, DICOM can only handle up to 9999 slices');
end

% Do a quick check to make sure that the number of slices matches the NIfTI header
if numberOfSlices ~= imageParamsNIfTI.ImageDim(3)
    error('Error: The number of slices in the image do not seem to match the number of slices specified in the NIfTI header!');
end

% Loop thru each slice in the image and write it to a DICOM file

k = waitbar(0,'Writing new DICOM series...');

imageDICOM.FileModDate=datestr(now);
imageDICOM.AcquisitionDate=datestr(now,'yyyymmdd');
imageDICOM.StudyDescription='Lead-DBS Plan';
imageDICOM.ProtocolName='Lead-DBS Plan';
imageDICOM.InstitutionName = 'Neuromodulation und Bewegungsstoerungen';
    
for j = 1:numberOfSlices
    
    %!! Extra: Get window tags off of each image
    fileName = get_filename(j);
    WindowGet = dicominfo([dicom_file]);
    imageDICOM.WindowWidth = WindowGet.WindowWidth;
    imageDICOM.WindowCenter = WindowGet.WindowCenter;

	% Update the DICOM header info
    imageDICOM.InstanceNumber = j;
    if ~isfield(imageDICOM,'SpacingBetweenSlices') % assume no gap.
        imageDICOM.SpacingBetweenSlices=zeros(3,1);
    end
    if j ~= 1
        imageDICOM.ImagePositionPatient(3) = imageDICOM.ImagePositionPatient(3) + imageDICOM.SpacingBetweenSlices(1);
        imageDICOM.SliceLocation           = imageDICOM.SliceLocation           + imageDICOM.SpacingBetweenSlices(1);
    end
    % Grab the current slice of the image
    outputImageSlice = outputImageMatrix(:,:,j);
    imageDICOM.LargestImagePixelValue  = max(max(outputImageSlice));
    imageDICOM.SmallestImagePixelValue = min(min(outputImageSlice));
    
    % load_nii pulls in the outputImageMatrix in RAS orientation.  Need to
    % flip the i and j axes to match the DICOM LPS orientation.
    % Also, the dicomwrite function in MATLAB seems to flip-flop the i/j
    % axes, so we need to do this also.
    outputReorientImageSlice = reorient_LPS(outputImageSlice);  %Need a new var name in caes the image size is not square
    
    % Finally, write to the DICOM
    fileName = get_filename(j);  % Generate the dicom filename
    dicomwrite(outputReorientImageSlice, fullfile(outputDirectory,fileName), imageDICOM, 'WritePrivate', true);
    % image(outputImageSlice, 'CDataMap', 'Scaled');
    
    waitbar(j/numberOfSlices);
end

if exist('k'); delete(k); end;  %#ok<EXIST> % Remove the waitbar
disp('Info: Output complete.');
return;



% <<Functions Below>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% -------------------------------------------------------------------------
% Get the image size and position parameters from the DICOM header
% Uses try/catch in case the headers are not defined in the DICOM file
function params = getDicomParams(imageDICOM, params)

try
    params.Orientation = imageDICOM.ImageOrientationPatient';
catch
    params.Orientation = [0 0 0 0 0 0];
    disp('Warning: DICOM: Orientation is not defined.');
end

temp = [0 0 1]; % Create a temporary vector to read in data
try
    temp(1) = imageDICOM.Rows;
    temp(2) = imageDICOM.Columns;
catch
    temp(1:2) = [0 0];
    disp('Warning: DICOM: the file you specified dosent seem to have the size');
    disp('         of its image array defined.');
end

try
    temp(3) = imageDICOM.ImagesInAcquisition;
catch
    temp(3) = 1;
    disp('Warning: DICOM: the file you specified seems to only be a single slice,');
    disp('         not a series of images.');
end

params.ImageDim = temp;

% Image Position
temp = [0 0 0]; % Voxel size - reset temp vector
try
    temp(1:2) = imageDICOM.PixelSpacing(1:2);
catch
    temp(1:2) = [0 0];
    disp('Warning: DICOM: Voxel dimensions are not defined');
end

try
    temp(3) = imageDICOM.SliceThickness;
catch
    try
        temp(3) = imageDICOM.SpacingBetweenSlices;
    catch
        temp(3) = 1;
        disp('Warning: DICOM: Slice thickness not defined.');
    end
end

params.VoxelDim = temp;

% Origin
try
    params.Origin = imageDICOM.ImagePositionPatient';   % Note that when we keep track of the
catch                                                   % Origin, we will always refer to the DICOM origin in LPS
    params.Origin = [0 0 0];
    disp('Warning: DICOM: Image origin not defined.');
end

return;

% -------------------------------------------------------------------------
% Get the image size and positon parameters from the NIfTI header
% We want everything represented in DICOM's LPS orientation, so some
% conversions will be done here
function params = getNiftiParams(imageNIfTI, params)

params.Orientation = [1 0 0 0 1 0];   % For now, we will always convert NIfTI into LPS
                                      % orientation, but in the future this setting may be changed
                                      % to allow any orientation to be specified

% Grab the image and voxel size from the NIfTI header
params.ImageDim = imageNIfTI.dim;
params.VoxelDim = imageNIfTI.voxsize;

% NIfTI uses the location of the upper-left (superior-left) corner as the origin  
% DICOM is in LPS orientation, so the x and y axes are flipped, so we need to 
% measure the position of the lower-right corner instead, and flip the sign

origin=imageNIfTI.mat*[0;0;0;1];

origin_vox = origin(1:3)';  % Get the location of the origin voxel
lowerright_vox = [0 0 1];
% Use ImageDim to get the position of LowerRight voxel
% Measure the length from the origin to the LowerRight voxel, then multiply
% by the voxel size to get the real-world coordinates:
lowerright_vox(1:2) = params.ImageDim(1:2); 
params.Origin = (lowerright_vox - origin_vox) .* params.VoxelDim;
% Flip the sign of x&y
params.Origin(1:2) = -params.Origin(1:2);


% See for more information on calculating orientation and origin:
% http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/quatern.html
return;

% -------------------------------------------------------------------------
function [] = show_userinfo(imageDICOM, imageParamsDICOM, imageNIfTI, imageParamsNIfTI)
disp('-------------------------------------------------------------------');
disp('DICOM Header Information:');
% Use a temp variable patientName to grab the patient's name:
try
    patientName = imageDICOM.PatientName.FamilyName;
catch
    patientName = 'Undefined';
end
disp(['Patient Name: ' patientName]);
try
    protocolName = imageDICOM.ProtocolName;
catch
    protocolName = ' ';
end
disp(['Protocol Name:  ' protocolName]);
disp(['Orientation: [' num2str(imageParamsDICOM.Orientation) ']']);
% Check if orientation matches LPS 'radiological'
if isequal(imageParamsDICOM.Orientation, [1 0 0 0 1 0])
    disp('Orientation: LPS, Axial Radiological Standard');
else
    disp('Orientation: is NOT LPS, the image in this DICOM may not be the same orientation as your new overlay image.');
end
disp(['Image Size: ' num2str(imageParamsDICOM.ImageDim)]);
disp(['Voxel Size: ' num2str(imageParamsDICOM.VoxelDim)]);
disp(['Origin:     ' num2str(imageParamsDICOM.Origin)]);
disp('------------------------------------------------------------------------');
disp('NIfTI Header Information:');
disp(['Image Description: LeadDBSExport']);
% Check the qform & sform codes for the type of coordinate system used:

% from http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html

%  /* [qs]form_code value: / / x,y,z coordinate system refers to: /
% /-----------------------/ /---------------------------------------*/
% #define NIFTI_XFORM_UNKNOWN 0       /! Arbitrary coordinates (Method 1). /
% #define NIFTI_XFORM_SCANNER_ANAT 1  /! Scanner-based anatomicalcoordinates /
% #define NIFTI_XFORM_ALIGNED_ANAT 2  /! Coordinates aligned to another file's, or to anatomical "truth". /
% #define NIFTI_XFORM_TALAIRACH 3     /! Coordinates aligned to Talairach-Tournoux Atlas; (0,0,0)=AC, etc. /
% #define NIFTI_XFORM_MNI_152 4       /! MNI 152 normalized coordinates. /

% If sform is not zero, xform_nii will reorient using the sform matrix:
% if imageNIfTI.original.hdr.hist.sform_code > 0
%     switch imageNIfTI.original.hdr.hist.sform_code
%         case 1
%             disp('Coordinate System: Scanner-based  [OK]');
%         case 2
%             disp('Coordinate System: Aligned to another file');
%             disp('WARNING: The coordinate system may not match your original scanner image coordinates!');
%         case 3
%             disp('Coordinate System: Talairach-Tournoux Atlas WARNING');
%             disp('WARNING: The coordinate system may not match your original scanner image coordinates!');
%         case 4
%             disp('Coordinate System: MNI 152 Normalized Coords Q');
%             disp('WARNING: The coordinate system may not match your original scanner image coordinates!');
%     end
% % Otherwise xform_nii will reorient using the qform matrix:
% elseif imageNIfTI.original.hdr.hist.sform_code > 0
%     switch imageNIfTI.original.hdr.hist.sform_code
%         case 1
%             disp('Coordinate System: Scanner-based');
%         case 2
%             disp('Coordinate System: Aligned to another file');
%             disp('WARNING: The coordinate system may not match your original scanner image coordinates!');
%         case 3
%             disp('Coordinate System: Talairach-Tournoux Atlas');
%             disp('WARNING: The coordinate system may not match your original scanner image coordinates!');
%         case 4
%             disp('Coordinate System: MNI 152 Normalized Coords');
%             disp('WARNING: The coordinate system may not match your original scanner image coordinates!');
%     end
% else
%     disp('Coordinate System: Unknown/Arbitrary');
%     disp('WARNING: The coordinate system may not match your original scanner image coordinates!');
% end

disp(['Orientation: [' num2str(imageParamsNIfTI.Orientation) ']']);
disp(['Orientation: LPS, Axial Radiological Standard']);   %#ok<NBRAK> % Always converted to LPS for now, in the future this may be changed to support other orientations 
disp(['Image Size: ' num2str(imageParamsNIfTI.ImageDim)]);
disp(['Voxel Size: ' num2str(imageParamsNIfTI.VoxelDim)]);
disp(['Origin:     ' num2str(imageParamsNIfTI.Origin)]);
disp('------------------------------------------------------------------------');

return;

% -------------------------------------------------------------------------
% Determine if the number is an integer.  Couldn't find this fucntion
% built in to MATLAB, but if you know of one please let me know!
function truefalse = check_integer(number)

truefalse = round(number) == number;
return;

% -------------------------------------------------------------------------
function outputImage = reorient_LPS(inputImage)

% outputImage = inputImage;  % Preallocate the outputImage array
% DONT DO THIS -- THIS IS BAD IF THE IMAGE ISNT SQUARE  (V1.1.1  Aug.24)

imageSize = size(inputImage);

for i = 1:imageSize(1)
    for j = 1:imageSize(2)
        I = imageSize(1) + 1 - i;
        J = imageSize(2) + 1 - j;
        % Flip the i and j indicies to match the way dicomwrite works in MATLAB
        outputImage(j,i) = inputImage(I,J);
    end
end

return;

% -------------------------------------------------------------------------
function fileName = get_filename(j)

if j < 10
    fileName = ['I000' num2str(j) '.dcm'];
elseif j < 100
    fileName = ['I00'  num2str(j) '.dcm'];
elseif j < 1000
    fileName = ['I0'   num2str(j) '.dcm'];
elseif j < 10000
    fileName = ['I'    num2str(j) '.dcm'];
else
    error('Error: Too many slices in the image, DICOM can only handle up to 9999 slices');
end
