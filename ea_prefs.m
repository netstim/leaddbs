function prefs=ea_prefs(patientname)

% determine preferences here. For filenames, the variable 'patientname' can
% be used in string-handling. This variable will be a string with the same name as the patient
% folder.

% load loaded prefs (-> prefs.lp)
try
load([fileparts(which('lead')),filesep,'ea_prefs']);
prefs.lp=lp;
end

%% general file handling:
prefs.prenii_unnormalized=['pre_tra.nii']; % not needed if schoenecker normalization is used.
prefs.tranii_unnormalized=['tra.nii'];
prefs.sagnii_unnormalized=['sag.nii'];
prefs.cornii_unnormalized=['cor.nii'];
prefs.ctnii_unnormalized=['fusion.nii'];
prefs.rawctnii_unnormalized=['CT.nii'];
prefs.ctnii_coregistered=['rCT.nii'];

prefs.patientdir=[patientname];
prefs.prenii=['lpre.nii'];
prefs.tranii=['ltra.nii'];
prefs.cornii=['lcor.nii'];
prefs.sagnii=['lsag.nii'];
prefs.ctnii=['lfusion.nii'];

prefs.gprenii=['glpre.nii'];
prefs.gtranii=['gltra.nii'];
prefs.gcornii=['glcor.nii'];
prefs.gsagnii=['glsag.nii'];
prefs.gctnii=['glfusion.nii'];



prefs.normmatrix=['lmat.txt'];


%% 2D-Export
prefs.d2.useprepost='post'; % can be 'post' or 'pre' to set the backdrop.

%% DICOM-Handling:

prefs.dicom.dicomfiles=0; % 1: delete DICOMs after conversion, 2: move DICOMs to pt/DICOM folder after conversion. 0: leave DICOMs where they were (not recommended: DICOMs will then be always re-imported from the import folder).

%% Normalization:
prefs.normalize.coreg='auto'; % set to 'manual' to include manual coregistration check steps.

%% DTI-files:

prefs.b0=['b0image.nii'];
prefs.FTR_unnormalized=['FTR.mat'];
prefs.FTR_normalized=['wFTR.mat'];
prefs.DTD=['DTD.mat'];
prefs.HARDI='HARDI.mat';
prefs.dti='dti.nii';
prefs.bval='dti.bval';
prefs.bvec='dti.bvec';
prefs.sampledtidicom='sample_dti_dicom.dcm'; % sample DICOM file of DTI image (used for trackvis export).

%% volumes:
prefs.hullmethod=2; % set 2 to use isosurface, 1 for concavehull, 0 for convexhull.
prefs.hullsmooth=5; % set to smooth hulldata. Only applies if isosurface is used. Only odd numbers allowed. Set to 0 if you don't want to smooth.
prefs.hullsimplify='auto'; % 0.1 would reduce hulldata to 10%. set to simplify hulldata. Set to 1 to not simplify. Only applies if isosurface is used.

%% labels:
prefs.lhullmethod=2; % set 2 to use isosurface, 1 for concavehull, 0 for convexhull.
prefs.lhullsmooth=3; % set to smooth hulldata. Only applies if isosurface is used. Only odd numbers allowed. Set to 0 if you don't want to smooth.
prefs.lhullsimplify='auto'; % 0.1 would reduce hulldata to 10% ? set to simplify hulldata. Set to 1 to not simplify. Only applies if isosurface is used.


%% fibers:

prefs.addfibers={}; % additional fibers to show.


%% lead server:

prefs.ls.autosave=0;
prefs.firstrun='off'; 
prefs.ls.dir='/Volumes/DBS_Lokas/bg_rest/lead_server/';

%% video export:
% for help see documentation of CaptureFigVid
% (http://www.mathworks.com/matlabcentral/fileexchange/41093-create-video-of-rotating-3d-plot/content/CaptureFigVid/CaptureFigVid.m)
prefs.video.path=[-20,10;-110,10;-190,80;-290,10;-380,10];
prefs.video.opts.FrameRate=24;
prefs.video.opts.Duration=20;
prefs.video.opts.Periodic=true;
