function prefs=ea_prefs(patientname)

% determine preferences here:

%% general file handling:
prefs.prenii_unnormalized=['pre_tra.nii']; % not needed if schönecker normalization is used.
prefs.tranii_unnormalized=['tra.nii'];
prefs.sagnii_unnormalized=['sag.nii'];
prefs.cornii_unnormalized=['cor.nii'];
prefs.ctnii_unnormalized=['fusion.nii'];
prefs.rawctnii_unnormalized=['CT.nii'];


prefs.patientdir=[patientname];
prefs.tranii=['ltra.nii'];
prefs.cornii=['lcor.nii'];
prefs.sagnii=['lsag.nii'];
prefs.ctnii=['lfusion.nii'];

prefs.gtranii=['gltra.nii'];
prefs.gcornii=['glcor.nii'];
prefs.gsagnii=['glsag.nii'];
prefs.gctnii=['glfusion.nii'];



prefs.normmatrix=['lmat.txt'];

%% DTI-files:

prefs.b0=['b0image.nii'];
prefs.FTR_unnormalized=['FTR.mat'];
prefs.FTR_normalized=['wFTR.mat'];

%% volumes:
prefs.hullmethod=2; % set 2 to use isosurface, 1 for concavehull, 0 for convexhull.
prefs.hullsmooth=5; % set to smooth hulldata. Only applies if isosurface is used. Only odd numbers allowed. Set to 0 if you don't want to smooth.
prefs.hullsimplify=0.15; %'auto'; %'auto'; %'auto'; %'auto'; %'auto'; %'auto'; % 0.1 would reduce hulldata to 10% ? set to simplify hulldata. Set to 1 to not simplify. Only applies if isosurface is used.

%% labels:
prefs.lhullmethod=2; % set 2 to use isosurface, 1 for concavehull, 0 for convexhull.
prefs.lhullsmooth=3; % set to smooth hulldata. Only applies if isosurface is used. Only odd numbers allowed. Set to 0 if you don't want to smooth.
prefs.lhullsimplify='auto'; % 0.1 would reduce hulldata to 10% ? set to simplify hulldata. Set to 1 to not simplify. Only applies if isosurface is used.


%% fibers:

prefs.addfibers={}; % additional fibers to show.


%% lead server:

prefs.ls.dir=[fileparts(which('lead')),filesep,'..',filesep,'lead_server',filesep];
prefs.ls.autosave=0;
prefs.firstrun='off'; 
