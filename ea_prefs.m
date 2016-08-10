function prefs=ea_prefs(patientname)

% determine preferences here. For filenames, the variable 'patientname' can
% be used in string-handling. This variable will be a string with the same name as the patient
% folder.

% load loaded prefs (-> prefs.lp)
try
    load([ea_getearoot,'ea_prefs']);
    prefs.lp=lp;
end

%% general settings:
prefs.pp.do=0; % use parallel processing if available.
prefs.pp.csize=4; % specify how many clusters to recruit.
prefs.pp.profile='local'; % specify parallel processing profile.

%% general file handling:
prefs.prenii_unnormalized='anat.nii'; % (preoperative) anatomical MR image
prefs.prenii_unnormalized_t1='anat_t1.nii'; % (preoperative) anatomical MR image
prefs.prenii_unnormalized_pd='anat_pd.nii'; % (preoperative) anatomical MR image
prefs.rawpreniis={prefs.prenii_unnormalized,prefs.prenii_unnormalized_t1,prefs.prenii_unnormalized_pd};

prefs.tranii_unnormalized='postop_tra.nii'; % postoperative primary anatomical MR image
prefs.sagnii_unnormalized='postop_sag.nii'; % postoperative optional anatomical MR image (saggital planes), for visualization only.
prefs.cornii_unnormalized='postop_cor.nii'; % postoperative optional anatomical MR image (coronar planes), for visualization only.
prefs.rawctnii_unnormalized='postop_ct.nii';
prefs.ctnii_coregistered='rpostop_ct.nii';

prefs.patientdir=patientname;
prefs.prenii='lanat.nii';
prefs.prenii_t1='lanat_t1.nii';
prefs.prenii_pd='lanat_pd.nii';
prefs.tranii='lpostop_tra.nii';
prefs.cornii='lpostop_cor.nii';
prefs.sagnii='lpostop_sag.nii';
prefs.ctnii='lpostop_ct.nii';

prefs.gprenii='glanat.nii';
prefs.gprenii_t1='glanat_t1.nii';
prefs.gprenii_pd='glanat_pd.nii';
prefs.gtranii='glpostop_tra.nii';
prefs.gcornii='glpostop_cor.nii';
prefs.gsagnii='glpostop_sag.nii';
prefs.gctnii='glpostop_ct.nii';





%% connectome files:
prefs.rest_prefix='res*.nii'; % raw resting state fMRI data search string
prefs.rest='rest.nii'; % default for dcm2nii export.
prefs.rest_default='rest.nii'; % default for dcm2nii export.

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


%% 2D-Export:
prefs.d2.useprepost='post'; % can be 'post' or 'pre' to set the backdrop.
prefs.d2.groupcolors='lead'; % can be 'maxdist' to use ea_distinguishable_colors by Timothy E. Holy, 'lead' to use a handpicked color set inspired by colorblender.com or 'lines' to use the matlab lines colormap (supports only seven colors).
prefs.d2.isovolsmoothed='s'; % set to 's' if you want to use the smoothed image and '' if you want to display the raw, unsmoothed version of the isovolume.
prefs.d2.isovolcolormap='jet'; % color map to use for plotting of isovolume heatmap.
prefs.d2.isovolsepcomb='combined'; % set to 'combined' to use the lr-combined isovolume and 'lr' to use the separate isovolumes.

%% 3D-Visualization:
prefs.d3.fiberstyle='tube'; % set to 'line' to show thin fibers
prefs.d3.fiberdiameter=0.1; % diameter of fibers ? only works in all ML versions with style==tube
prefs.d3.maxfibers=200; % set to inf to show all fibers (but this could lead to crashes).

%% DICOM-Handling:

prefs.dicom.dicomfiles=1; % 1: delete DICOMs after conversion, 0: Leave DICOMs at pt/DICOM folder after conversion.

%% Normalization:
prefs.normalize.coreg='auto'; % set to 'manual' to include manual coregistration check steps.
prefs.normalize.default='ea_normalize_ants';
prefs.normalize.inverse.customtpm=0; % set to 1 if custom TPM shall be built for inverse warpings.

%% Coregistration (CT/MR):
prefs.ctcoreg.default='ea_coregctmri_ants';

%% Coregistration (MR/MR):
prefs.mrcoreg.default='spm'; % set to 'spm' or 'ants'

%% volumes:
prefs.hullmethod=2; % set 2 to use isosurface, 1 for concavehull, 0 for convexhull.
prefs.hullsmooth=5; % set to smooth hulldata. Only applies if isosurface is used. Only odd numbers allowed. Set to 0 if you don't want to smooth.
prefs.hullsimplify=0.5; % 0.1 would reduce hulldata to 10%. set to simplify hulldata. Set to 1 to not simplify. Only applies if isosurface is used.

%% labels:
prefs.lhullmethod=2; % set 2 to use isosurface, 1 for concavehull, 0 for convexhull.
prefs.lhullsmooth=7; % set to smooth hulldata. Only applies if isosurface is used. Only odd numbers allowed. Set to 0 if you don't want to smooth.
prefs.lhullsimplify='auto'; % 0.1 would reduce hulldata to 10% ? set to simplify hulldata. Set to 1 to not simplify. Only applies if isosurface is used.


%% DICOM files:
prefs.dicom.dicomfiles=0; % 1: delete DICOMs after conversion, 2: move DICOMs to pt/DICOM folder after conversion. 0: leave DICOMs where they were (not recommended: DICOMs will then be always re-imported from the import folder).


%% fibers:

prefs.addfibers={}; % additional fibers to show.


%% lead server:

prefs.ls.autosave=0;

%% xelatex executable path:
prefs.ltx.pdfconverter='/Library/TeX/texbin/xelatex';

%% environment
prefs.env.dev=0;
prefs.firstrun='off';

prefs.ls.dir='/PA/Neuro/_projects/lead/lead_server/';
prefs.ixi.dir='/Volumes/TimeMachine/IXI/DB/';