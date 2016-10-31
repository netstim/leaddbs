% NAME: mrAlign (version 4.3)
% AUTHOR: DJH
% DATE: 6/2004
% PURPOSE:
%	Aligning/registering inplane and volume anatomies.
% HISTORY:
%	Dozens of people contributed to the original mrAlign. It was written
%   using scripts in matlab 4 without any data structures. mrAlign4  
%   cleaned up the code. Version 4.1 adds a GUI. Version 4.2 ported to
%   matlab 7 and adds a number of additional features.
%        $Id$


%%%%%%%%%%%%%%%%%%%
% Global Variables %
%%%%%%%%%%%%%%%%%%%

global ALIGN
global mrDEFAULTS

% Load .mrDefaults
mrDEFAULTS = loadMrDefaults;

% Check Matlab version number
[mlrVersion, expectedMatlabVersion] = mrLoadRetVersion;
version = ver('Matlab');
matlabVersion = str2num(version.Version(1:3));
if ~ismember(matlabVersion, expectedMatlabVersion);
    mrWarnDlg(['mrAlign is intended for Matlab ',num2str(expectedMatlabVersion),...
        '. You are running Matlab ',version.Version]);
end

ALIGN.version = mlrVersion;

ALIGN.volumePath = [];           % path string to volume anatomy file
ALIGN.volume = [];               % volume matrix
ALIGN.volumeHdr = [];            % loaded from nifti header
ALIGN.volSize = [];              % size of volume matrix
ALIGN.volumeVoxelSize = [];      % 3x1 vector specifying voxel size
ALIGN.volumePermutation = eye(3);  % extracted from volume header to determine slice orientation
ALIGN.volumeClip = [];           % for display

ALIGN.inplanePath = [];          % path string to inplane anatomy file
ALIGN.inplanes = [];             % inplane matrix
ALIGN.inplaneHdr = [];           % loaded from nifti header
ALIGN.inplaneSize = [];          % size of inplane matrix
ALIGN.inplaneVoxelSize = [];     % 3x1 vector specifying voxel size
ALIGN.inplanesClip = [];         % for display

ALIGN.crop = [];                 % specifies 2D crop region used for contrast correction
ALIGN.guiXform = [];             % 4x4 transform matrix read from rot and trans GUI
ALIGN.xform = [];                % 4x4 transform matrix from inplane pixels -> volume pixels. this is the current transformation
ALIGN.xformICCorrection = [];    % 4x4 transform matrix from inplane pixels -> volume pixels. this is the transformation after the last intensity/contrast correction was computed
ALIGN.oldXform = [];             % 4x4 transform matrix from inplane pixels -> volume pixels. this is the transformation before computing the last alignment
ALIGN.oldGuiRotation = [];       % 3-tuple to keep old manual rotation in memory after running alignment
ALIGN.oldGuiTranlastion = [];    % 3-tuple to keep old manual translation in memory after running alignment
ALIGN.NIter = 200;               % Number of iterations in auto alignment
ALIGN.sliceOrientation = 1;      % 1=axial, 2=coronal, 3=sagittal
ALIGN.coords = [1,1,1];          % selected voxel
ALIGN.currentlyComputingAlignment = 0;% an alignment is currently being computed
ALIGN.stopComputingAlignment = 0;     % computing alignment is being cancelled

ALIGN.robustAlignment = 1;       % specifies if robustMest is used to optimize the alignment 
ALIGN.rigidBodyAlignment = 1;    % specifies if alignement is restricted to 6 or 7 degrees of freedom
ALIGN.reverseContrastAlignment = 0; % specifies if contrast is reversed for alignment optimization 
ALIGN.ignoreZeroVoxels = 1;      % specifies if zero voxels should be ignored when optimizing the alignment 
ALIGN.displayAlignmentSteps = 1; % specifies if alignment optimization steps are displayed

ALIGN.minCoarseVoxelSize = 3;    % minimum size of inplanes/volume voxels to compute coarse alignment
ALIGN.minFineVoxelSize = 0.5;      % minimum size of inplanes voxels to compute fine alignment
ALIGN.minEPICoarseVoxelSize = 2; % minimum size of inplanes/volume voxels to compute coarse EPI alignment
ALIGN.minEPIFineVoxelSize = .5;  % minimum size of inplanes voxels to compute fine EPI alignment

ALIGN.correctedInplanes = [];    % to store the intensity/contrast corrected inplanes
ALIGN.correctedVolume = [];      % to store the intensity/contrast corrected volume
ALIGN.reversedContrast = 0;      % specifies if current corrected volumes include reversed contrast
ALIGN.xformIsRigidBody = 1;      % specifies if current xform is rigid body or not


ALIGN.baseCmap = gray(256);
ALIGN.overlayCmap = hot(272);
ALIGN.overlayCmap = ALIGN.overlayCmap(17:end,:); %let's remove a bit of black to tell apart zeros from source and destination

disp(['mrAlign (version ',num2str(ALIGN.version),')']);

mrAlignGUI


clear  mlrVersion expectedMatlabVersion version matlabVersion

return
