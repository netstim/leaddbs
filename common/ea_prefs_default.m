function prefs = ea_prefs_default(patientname)

% determine preferences here. For filenames, the variable 'patientname' can
% be used in string-handling. This variable will be a string with the same name as the patient
% folder.

% load loaded prefs (-> prefs.lp)
try
    load([ea_getearoot,'ea_prefs']);
    prefs.lp = lp;
end

%% adjust to user
prefs.dev.profile='user';

%% general settings:
prefs.pp.do=0; % use parallel processing if available.
prefs.pp.csize=4; % specify how many clusters to recruit.
prefs.pp.profile='local'; % specify parallel processing profile.

%% general file handling:
prefs.prenii_searchstring='anat_*.nii';
prefs.prenii_order={'t1','t2','pd'}; % assign order of anatomical images used (nonspecified images matching the searchstring will be added afterwards).
prefs.prenii_unnormalized='anat_t2.nii'; % default (still used for DICOM import)
prefs.prenii_unnormalized_t1='anat_t1.nii'; % default (still used for DICOM import)
prefs.prenii_unnormalized_pd='anat_pd.nii'; % default (still used for DICOM import)

prefs.tranii_unnormalized='postop_tra.nii';
prefs.sagnii_unnormalized='postop_sag.nii';
prefs.cornii_unnormalized='postop_cor.nii';
prefs.rawctnii_unnormalized='postop_ct.nii';
prefs.ctnii_coregistered='rpostop_ct.nii';
prefs.tp_ctnii_coregistered=['tp_',prefs.ctnii_coregistered];

prefs.preferMRCT = 1; % preference of MR or CT modality for post-op image: 1 for MR, 2 for CT.

prefs.patientdir=patientname;

prefs.gprenii='glanat.nii';
prefs.gtranii='glpostop_tra.nii';
prefs.gcornii='glpostop_cor.nii';
prefs.gsagnii='glpostop_sag.nii';
prefs.gctnii='glpostop_ct.nii';

prefs.tp_gctnii=['tp_',prefs.gctnii];

%% BIDS settings
prefs.bids_session_postop = 'ses-postDBS';
prefs.bids_session_preop = 'ses-preDBS';

%% Misc:
prefs.tonemap='heuristic'; % set to 'albada' to change to datadriven mode

%% connectome files:
prefs.rest_searchstring='rest*.nii'; % raw resting state fMRI data search string
prefs.rest='rest.nii'; % default for dcm2nii export.

%% connectome settings:
prefs.lc.struc.maxdist=2; % maximal distance to form a connection (between fiber terminals and voxel centers, in mm).
prefs.lc.struc.minlen=3; % minimum fiber length to consider to form connections (in segments).
prefs.lc.graphsurfc=[0.2081 0.1663 0.5292]; % default color for graph metric 3D-visualizations.
prefs.lc.matsurfc=[0.8 0.7 0.4]; % default color for matrix-level correlations 3D-visualizations.
prefs.lc.seedsurfc=[0.8 0.1 0.1]; % default color for seed of matrix-level correlations 3D-visualizations.
prefs.lc.func.regress_global=1;
prefs.lc.func.regress_wmcsf=1;
prefs.lc.func.bphighcutoff=0.08;
prefs.lc.func.bplowcutoff=0.009;
prefs.lc.datadir=[ea_getearoot,'connectomes',filesep];

%% connectome mapper settings:
prefs.lcm.vatseed='binary'; % set to 'efield_gauss' to use weighted seed of normalized E-field (or 'efield' to use weighted seed of unmodified e-field - not recommended).
prefs.lcm.chunk=10; % define how many fMRI seeds to handle in the same run. Can be 0 to handle all supplied. Depending on RAM available, 5-20 is a good option.
prefs.lcm.includesurf=1; % if surface definitions are available for connectomes, include those, too
prefs.lcm.struc.patienttracts.nativeseed = 0; % Use native space VAT for calculation when patient's fiber tracts is selected.

%% DTI-files:
prefs.b0='b0.nii';
prefs.fa='fa.nii';
prefs.fa2anat='fa2anat.nii';
prefs.FTR_unnormalized='FTR.mat';
prefs.FTR_normalized='wFTR.mat';
prefs.DTD='DTD.mat';
prefs.HARDI='HARDI.mat';
prefs.dti='dti.nii';
prefs.bval='dti.bval';
prefs.bvec='dti.bvec';
prefs.sampledtidicom='sample_dti_dicom.dcm'; % sample DICOM file of DTI image (used for trackvis export).

prefs.normmatrix='lmat.txt';

%% Normalization:
% prefs.normalize.coreg='auto'; % set to 'manual' to include manual coregistration check steps.
prefs.normalize.default='ea_normalize_ants';
prefs.normalize.inverse.warp='inverse'; % set to 'tpm' in case you wish to create a atlas-specific tpm to warp atlases, set to 'inverse' to apply the inverse transform of your normalization.
prefs.normalize.inverse.customtpm=0; % set to 1 if custom TPM shall be built for inverse warpings. Only applies if the above is set to 'tpm'.
prefs.normalize.createwarpgrids=0; % set to 1 to create grid files that show deformation fields in "Show Normalization" option.
prefs.normalize.fsl.warpres=8; % Defines the warp resolution in FSL warps.
prefs.normalize.spm.resolution=1; % Defines resolution in mm when using SPM normalization routines (New Segment, DARTEL, SHOOT).

%% Reconstruction
prefs.reco.mancoruse='postop'; % switch to 'rpostop' to use resliced CT.
prefs.reco.saveACPC=0; % also save reconstructions in AC/PC space
prefs.reco.saveimg=0; % save fiducial marker visualisation as image after "Refined TRAC/CORE"
prefs.reco.exportfiducials=0; % automatically export fiducials to a comma separated value file after "Refined TRAC/CORE". Set this to 'fcsv' for simple import into Slicer, otherwise set to 'csv' or 'txt' for import into other software.

%% Coregistration (CT/MR):
prefs.ctcoreg.default='ea_coregctmri_ants';

%% Coregistration (MR/MR):
prefs.mrcoreg.default='spm'; % set to 'spm' or 'ants'
prefs.mrcoreg.writeoutcoreg=0; % set default to 0 to prevent writing out coregistration transformation

%% Subcortical refine (Post to Pre):
prefs.scrf.tonemap='tp_'; % can set to '' if want to use non-tonemapped CTs for brainshift correction (default = 'tp_').

%% Atlas-Sets:
prefs.atlases.default='DISTAL Minimal (Ewert 2017)';

%% volumes:
prefs.hullmethod=2; % set 2 to use isosurface, 1 for concavehull, 0 for convexhull.
prefs.hullsmooth=5; % set to smooth hulldata. Only applies if isosurface is used. Only odd numbers allowed. Set to 0 if you don't want to smooth.
prefs.hullsimplify=0.3; % 0.1 would reduce hulldata to 10%. set to simplify hulldata. Set to 1 to not simplify. Only applies if isosurface is used.

%% labels:
prefs.lhullmethod=2; % set 2 to use isosurface, 1 for concavehull, 0 for convexhull.
prefs.lhullsmooth=3; % set to smooth hulldata. Only applies if isosurface is used. Only odd numbers allowed. Set to 0 if you don't want to smooth.
prefs.lhullsimplify='auto'; % 0.1 would reduce hulldata to 10% ? set to simplify hulldata. Set to 1 to not simplify. Only applies if isosurface is used.

%% 2D-Export
prefs.d2.useprepost='pre'; % can be 'post' or 'pre' to set the backdrop.
prefs.d2.groupcolors='lead'; % can be 'maxdist' to use ea_distinguishable_colors by Timothy E. Holy, 'lead' to use a handpicked color set inspired by colorblender.com or 'lines' to use the matlab lines colormap (supports only seven colors).
prefs.d2.isovolsmoothed='s'; % set to '' if you want to display the raw, unsmoothed version of the isovolume.
prefs.d2.isovolcolormap='jet'; % color map to use for plotting of isovolume heatmap.
prefs.d2.isovolsepcomb='combined'; % set to 'combined' to use the lr-combined isovolume and 'lr' to use the separate isovolumes.

%% 3D-Visualization:
prefs.d3.fiberstyle='tube'; % set to 'line' to show thin fibers
prefs.d3.fiberdiameter=0.1; % diameter of fibers ? only works in all ML versions with style==tube
prefs.d3.maxfibers=200; % set to inf to show all fibers (but this could lead to crashes).
prefs.d3.colorjitter=0; % set to 0 to show no color jitter at all.
prefs.d3.showdirarrows = 0;

%% Video export
prefs.video.path=[-90,10;-110,10;-180,80;-250,10;-360,10;-450,10];
prefs.video.opts.FrameRate=24;
prefs.video.opts.Duration=30;
prefs.video.opts.Periodic=true;

%% FEM-VAT settings:
prefs.vat.gm='mask'; % set to 'atlas' to use current atlas single structures, 'mask' to use 'gm_mask.nii', set to 'tpm' to use c1 portion of tpm.
prefs.vat.efieldmax=10000; % set upper limit for maximal values in efield.

%% MER-Visualization:
prefs.mer.rejwin = [1 60];
prefs.mer.offset = 2; % default distance between mer tracts is 2mm
prefs.mer.length = 24; % default mer length for visualization is 24mm
prefs.mer.markersize = 0.5; % default mer marker size 0.25mm
prefs.mer.defaulttract = 1; % default tract is Central(1). Set to 2=Anterior,3=Posterior,4=Lateral, or 5=Medial
prefs.mer.n_pnts = 50;
prefs.mer.tag.visible = 'off';
prefs.mer.step_size = [0.25 0.75 0.05];
prefs.mer.tract_info = struct(...
    'label', {'central', 'anterior', 'posterior', 'lateral', 'medial'},...
    'color', {  [0.5,0,0],...       Maroon
                [0.5,0.5,0],...     Olive
                [0,0.5,0],...       Green
                [0.5,0,0.5],...     Purple
                [0,0.5,0.5]},...    Teal; 0,0,0.5 Navy
    'position', {  [ 0,  0, 0],...
                    [ 0,  1, 0],...
                    [ 0, -1, 0],...
                    [ 1,  0, 0],...
                    [-1,  0, 0]});

%% Cortex-Visualization:
prefs.d3.cortexcolor=[0.65 0.65 0.65]; % default color is gray
prefs.d3.cortexalpha=0.5; % default alpha is 0.5
prefs.d3.cortex_defaultatlas='DKT'; % Currently supports 'DKT','DKT_aseg','a2009'

%% Freesurfer Preferences
prefs.d3.fs.dev=0;

%% DICOM files:
prefs.dicom.dicomfiles=0; % 1: delete DICOMs after conversion, 0: Leave DICOMs at pt/DICOM folder after conversion.
prefs.dicom.tool='dcm2niix'; % switch to 'dicm2nii' to use the Matlab based converter.

%% fibers:

prefs.addfibers={}; % additional fibers to show.

%% native-space:
prefs.native.warp='inverse'; % set to 'tpm' in case you wish to create a atlas-specific tpm to warp atlases, set to 'inverse' to apply the inverse transform of your normalization.

%% lead server:
prefs.ls.autosave=0;

%% environment
prefs.env.dev=0;
prefs.env.logtime=0;
prefs.env.campus='generic';
prefs.ixi.meanage=60; % mean age used if no patient/subject age is specified in folder.

%% xelatex executable path:
prefs.ltx.pdfconverter=''; % set path to xelatex here (for PDF export)

prefs.ls.dir=''; % set path to lead server here (for web export)
prefs.ixi.dir=''; % set path to ixi database here
prefs.ixi.meanage=60; % mean age used if no patient/subject age is specified in folder.


%% genetics
prefs.genetics.dbdir=[ea_space,'genetics',filesep];


%% platform specific

% Set to true this line if libstdc++.so.6 is needed.
% However it is preferrable to fix it at system level (e.g. using package build-essentials). 
% Additionally, install the matlab-support package and choose to use the system libraries for gcc.
% If set to true it will add this path: addpath(fullfile(earoot,'ext_libs\support\unix'));
prefs.platform.unix.load_shipped_libstdcpp6=false; 
