% configure file for operations on streamline trackings results
%
% for BATCH EDITOR system (Volkmar Glauche)
%
% File created by Susanne Schnell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ftrop = dti_cfg_ftrop

% ---------------------------------------------------------------------
% filename1 name of ftrstruct
% ---------------------------------------------------------------------
filename1         = cfg_files;
filename1.tag     = 'filename1';
filename1.name    = 'Select the ftrstruct';
filename1.help    = {'Select the ftrstruct (usually something like blah_FTR.mat'};
filename1.filter  = 'mat';
filename1.ufilter = '.*';
filename1.num     = [1 1];

% ---------------------------------------------------------------------
% filename2 name of dtdstruct
% ---------------------------------------------------------------------
filename2         = cfg_files;
filename2.tag     = 'filename2';
filename2.name    = 'Select the dtdstruct';
filename2.help    = {'Select the dtdstruct (usually something like "blah_DTD.mat"'};
filename2.filter  = '_DTD.mat';
filename2.ufilter = '.*';
filename2.num     = [1 1];

% ---------------------------------------------------------------------
% filename2 name of dtdstruct
% ---------------------------------------------------------------------
filename3         = cfg_files;
filename3.tag     = 'filename3';
filename3.name    = 'Select reference frame';
filename3.help    = {'A mrstruct, or a dtdstruct, or a nifti determining the volume of interest.'};
filename3.filter  = '.mat|.nii';
filename3.ufilter = '.*';
filename3.num     = [1 1];

% ---------------------------------------------------------------------
% filename2 name of dtdstruct
% ---------------------------------------------------------------------
filename4         = cfg_files;
filename4.tag     = 'filename4';
filename4.name    = 'Select Contrasts';
filename4.help    = {'A mrstruct, or a dtdstruct, or a nifti used for tractometry'};
filename4.filter  = '.mat|.nii';
filename4.ufilter = '.*';
filename4.num     = [1 inf];


% ---------------------------------------------------------------------
% filename2 name of dtdstruct
% ---------------------------------------------------------------------
filename5         = cfg_files;
filename5.tag     = 'filename5';
filename5.name    = 'Select ROIs';
filename5.help    = {'A maskstruct, mrsrtuct or a nifti containg ROIs'};
filename5.filter  = '.mat|.nii';
filename5.ufilter = '.*';
filename5.num     = [1 inf];



% ---------------------------------------------------------------------
% filedef name of deformation field
% ---------------------------------------------------------------------
filedef         = cfg_files;
filedef.tag     = 'filedef';
filedef.name    = 'Select the deformation field';
filedef.help    = {'Select the deformation field (usually created by "New Segment", something like "yblah.nii"'};
filedef.filter  = '.nii';
filedef.ufilter = '.*';
filedef.num     = [1 1];

% ---------------------------------------------------------------------
% filedefi name of inverse deformation field
% ---------------------------------------------------------------------
filedefi         = cfg_files;
filedefi.tag     = 'filedefi';
filedefi.name    = 'Select the inverse deformation field';
filedefi.help    = {'Select the inverse deformation field (usually created by "New Segment", something like "iyblah.nii"'};
filedefi.filter  = '.nii';
filedefi.ufilter = '.*';
filedefi.num     = [1 1];


% ---------------------------------------------------------------------
% maskname, select the correct mask by name
% ---------------------------------------------------------------------
maskname         = cfg_entry;
maskname.tag     = 'maskname';
maskname.name    = 'Type in the name of the mask';
maskname.help    = {'Type in the name of the mask with which the fibers will be selected.'};
maskname.strtype = 's';
maskname.num     = [1 Inf];

% ---------------------------------------------------------------------
% masknumber, select the correct mask by a number
% ---------------------------------------------------------------------
masknumber        = cfg_entry;
masknumber.tag     = 'masknumber';
masknumber.name    = 'Mask Number';
masknumber.help    = {'Select the mask by a number. You need to remember the mask position by yourself.'};
masknumber.strtype = 'e';
masknumber.num     = [1 1];

% ---------------------------------------------------------------------
% mask, select the correct mask by a number
% ---------------------------------------------------------------------
mask         = cfg_choice;
mask.tag     = 'mask';
mask.name    = 'Mask Number or Mask Name';
mask.help    = {'Select the mask by a number or name. Please take care to spell the name correctly.'...
                'In case of mask number, you need to remember the mask position by yourself.'};
mask.values  = {masknumber maskname};

% ---------------------------------------------------------------------
% roiname
% ---------------------------------------------------------------------
roiname         = cfg_files;
roiname.tag     = 'roiname';
roiname.name    = 'Load ROI';
roiname.help    = {'Select a maskstruct (ROI).'};
roiname.filter  = 'mat';
roiname.ufilter = '.*';
roiname.num     = [1 1];

% ---------------------------------------------------------------------
% roidef 
% ---------------------------------------------------------------------
roidef         = cfg_branch;
roidef.tag     = 'roidef';
roidef.name    = 'Definition of ROI';
roidef.help    = {'Choose the maskstruct-file (ROI) and select the mask.'};
roidef.val     = {roiname mask};

% ---------------------------------------------------------------------
% newfilename name of resulting file
% ---------------------------------------------------------------------
newfilename         = cfg_entry;
newfilename.tag     = 'newfilename';
newfilename.name    = 'New File Name';
newfilename.val     = {'.mat'};
newfilename.help    = {'Type in the name of the new file, if you leave this empty a default file name is generated with the endings "_selFTR.mat", "_defFTR.mat", "_addMAP.mat" or "_multMAP.mat"'};
newfilename.strtype = 's';
newfilename.num     = [1 Inf];


% ---------------------------------------------------------------------
% newfilename name of resulting file
% ---------------------------------------------------------------------
newfilename2         = cfg_entry;
newfilename2.tag     = 'newfilename2';
newfilename2.name    = 'Output File Name';
newfilename2.help    = {'Type in the name of the file, where statistics is saved'};
newfilename2.strtype = 's';
newfilename2.num     = [1 Inf];


% ---------------------------------------------------------------------
% fibername name of new fiber subset
% ---------------------------------------------------------------------
fibername         = cfg_entry;
fibername.tag     = 'fibername';
fibername.name    = 'Name of fiber subset';
fibername.val     = {'fiber_001'};
fibername.help    = {'Type in the name of the new fiber subset'};
fibername.strtype = 's';
fibername.num     = [1 Inf];

% ---------------------------------------------------------------------
% thresh, select the threshhold for creating a mask
% ---------------------------------------------------------------------
thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Threshold for creating a mask';
thresh.val     = {0};
thresh.help    = {'Threshold for creating a mask of voxels visited by at least this number of fibretracks. Default = 0'};
thresh.strtype = 'e';
thresh.num     = [1 1];

% ---------------------------------------------------------------------
% fibernumber, select the correct fibersubset by a number
% ---------------------------------------------------------------------
fibernumber         = cfg_entry;
fibernumber.tag     = 'fibernumber';
fibernumber.name    = 'Fibersubset Number';
fibernumber.val     = {1};
fibernumber.help    = {'Select the fiber subset by a number. You need to remember the position by yourself.'};
fibernumber.strtype = 'e';
fibernumber.num     = [1 1];


% ---------------------------------------------------------------------
% oversampling
% ---------------------------------------------------------------------
osamp        = cfg_entry;
osamp.tag     = 'oversampling';
osamp.name    = 'Oversampling Factor';
osamp.val     = {1};
osamp.help    = {'Select the oversampling factor. The voxelsize of the reference frame are divied by this factor.'};
osamp.strtype = 'e';
osamp.num     = [1 1];

% ---------------------------------------------------------------------
% fibersubset, choose whether to select the fiber subset number or name
% ---------------------------------------------------------------------
fibersubset         = cfg_choice;
fibersubset.tag     = 'fibersubset';
fibersubset.name    = 'Fibersubset Number or Name';
fibersubset.help    = {'Select the fiber subset by a number or by name. You need to remember the correct spelling or the position by yourself.'};
fibersubset.values  = {fibernumber fibername};

% ---------------------------------------------------------------------
% visitMask input only the ftrstruct
% ---------------------------------------------------------------------
endpointMask         = cfg_exbranch;
endpointMask.tag     = 'endpointMask';
endpointMask.name    = 'Create end point mask';
endpointMask.help    = {'Determines the endpoints of the fibertracks in the specified subset and creates a mask of it'
    'Select the ftrstruct and the fiber subset.'};
endpointMask.val     = {filename1 fibersubset filename2 newfilename};
endpointMask.prog    = @(job)dti_ftrop_ui('endpointMask',job);
endpointMask.vout    = @vout;

% ---------------------------------------------------------------------
% visitMask input the ftrstruct and a threshold
% ---------------------------------------------------------------------
visitMask         = cfg_exbranch;
visitMask.tag     = 'visitMask';
visitMask.name    = 'Create a Visit Mask';
visitMask.help    = {'Creates a mask of voxels visited by the fibertracks of the specified subset.'
    'A single fibertrack visiting a voxel twice, will increment the voxel also by two. Select the ftrstruct and the fiber subset.'
    'The resulting mask will be saved as maskstruct and can be loaded in "fiberviewer_tool" instead as ROI'
    'by using the -> open -> ROI dialogue'};
visitMask.val     = {filename1 fibersubset thresh filename2 newfilename};
visitMask.prog    = @(job)dti_ftrop_ui('visitMask',job);
visitMask.vout    = @vout;

% ---------------------------------------------------------------------
% visitMap input only the ftrstruct
% ---------------------------------------------------------------------
visitMap         = cfg_exbranch;
visitMap.tag     = 'visitMap';
visitMap.name    = 'Create Visit Map';
visitMapext.help    = {'Determines how often a voxel was visited by the fibertracks containing in the specified subset.'
    'A single fibertrack visiting a voxel twice, will increment the voxel also by two. Select the ftrstruct and the fiber subset.'
    'The resulting map will be saved as mrstruct and can be loaded in "fiberviewer_tool" instead of dtdstruct or in addition to'
    'dtdstruct buy using the -> open -> append volume dialogue'};
visitMap.val     = {filename1 fibersubset filename2 newfilename};
visitMap.prog    = @(job)dti_ftrop_ui('visitMap',job);
visitMap.vout    = @vout;

% ---------------------------------------------------------------------
% visitMap input only the ftrstruct
% ---------------------------------------------------------------------
visitMapext         = cfg_exbranch;
visitMapext.tag     = 'visitMapext';
visitMapext.name    = 'Create Visit Map (extended)';
visitMap.help    = {'Determines how often a voxel was visited by the fibertracks containing in the specified subset.'
    'The tracts are rendered into the volume by trilinear interpolation. Three files are created: a color-coded tract density, and'
    'a ordinaray gray-scale density, and endpoint densitity. Files may saved as mrstruct or nifits (extension of destination file determines type)'
    'The reference image (nifti or mrstruct) determines volume of interest and voxelsizes.'};
visitMapext.val     = {filename1 fibersubset filename3 osamp newfilename};
visitMapext.prog    = @(job)dti_ftrop_ui('visitMapext',job);
visitMapext.vout    = @voutFDext;


% ---------------------------------------------------------------------
% eliminatebyROI input of ftrstruct and ROI
% ---------------------------------------------------------------------
eliminatebyROI         = cfg_exbranch;
eliminatebyROI.tag     = 'eliminatebyROI';
eliminatebyROI.name    = 'Eliminate fibers by ROI';
eliminatebyROI.help    = {'Select the two maps you want to multiply. Multiplaction results depend on the probabilistic tracking method chosen. If the extended version was used you will have three maps: the merging fibre map, the connecting fibre map and all-in-one map.'};
eliminatebyROI.val     = {filename1 fibersubset roidef fibername newfilename};
eliminatebyROI.prog    = @(job)dti_ftrop_ui('eliminatebyROI',job);
eliminatebyROI.vout    = @vout;



seltype         = cfg_menu;
seltype.tag     = 'seltype';
seltype.name    = 'Selection Method';
seltype.help    = {'Choose the way the fiber is defined to belong to a certain ROI.'};
seltype.labels = {
                 'visiting the area'
                 'endpoint in Area'
}';
seltype.values = {0 1};
seltype.val = {0};




% ---------------------------------------------------------------------
% selectbyROI input of one ftrstruct, a maskstruct and the new filename
% ---------------------------------------------------------------------
selectbyROI         = cfg_exbranch;
selectbyROI.tag     = 'selectbyROI';
selectbyROI.name    = 'Select Fiber subsets by ROIs';
selectbyROI.help    = {'Select fiber subsets with help of ROIs (maskstructs).'};
selectbyROI.val     = {filename1 fibersubset roidef seltype fibername newfilename};
selectbyROI.prog    = @(job)dti_ftrop_ui('selectbyROI',job);
selectbyROI.vout    = @vout;




% ---------------------------------------------------------------------
% selectbyROI input of one ftrstruct, a maskstruct and the new filename
% ---------------------------------------------------------------------


overlapratio         = cfg_entry;
overlapratio.tag     = 'overlapratio';
overlapratio.name    = 'Overlap ratio';
overlapratio.val     = {0.1};
overlapratio.help    = {'fibers that are have at least this ratio of overlap with the ROI are selected'};
overlapratio.strtype = 'e';
overlapratio.num     = [1 1];



selectbyOverlapROI         = cfg_exbranch;
selectbyOverlapROI.tag     = 'selectbyOverlapROI';
selectbyOverlapROI.name    = 'Select Fiber subsets by overlap with ROIs';
selectbyOverlapROI.help    = {'Select fiber subsets with help of ROIs (maskstructs).'};
selectbyOverlapROI.val     = {filename1 fibersubset roidef overlapratio fibername newfilename};
selectbyOverlapROI.prog    = @(job)dti_ftrop_ui('selectbyOverlapROI',job);
selectbyOverlapROI.vout    = @vout;




% ---------------------------------------------------------------------
% fuzzy or not fuzzy ? 
% ---------------------------------------------------------------------

threshfuz         = cfg_entry;
threshfuz.tag     = 'threshfuz';
threshfuz.name    = 'Threshold for selection';
threshfuz.val     = {0.1};
threshfuz.help    = {'threshold'};
threshfuz.strtype = 'e';
threshfuz.num     = [1 1];

sigmafuz         = cfg_entry;
sigmafuz.tag     = 'sigmafuz';
sigmafuz.name    = 'Sigma of Gaussian Smooth';
sigmafuz.val     = {1};
sigmafuz.help    = {'threshold'};
sigmafuz.strtype = 'e';
sigmafuz.num     = [1 1];

fuzzy         = cfg_branch;
fuzzy.tag     = 'fuzzy';
fuzzy.name    = 'Select parameters for fuzzy selection';
fuzzy.val     = {threshfuz sigmafuz};
fuzzy.help    = {'Specify a output filename.'};

nofuzzy        = cfg_const;
nofuzzy.tag    = 'nofuzzy';
nofuzzy.name   = 'No fuzzy selection';
nofuzzy.val    = {true};
nofuzzy.help   = {'Select streamlines with hard ROI boundaries.'};

fuzzysel         = cfg_choice;
fuzzysel.tag     = 'fuzzysel';
fuzzysel.name    = 'Use fuzzy Selection';
fuzzysel.values  = {fuzzy nofuzzy};
fuzzysel.help    = {'Use fuzzy selection to get more stabe fiber counts'};


%%%%%%%%%% contrast for tractometry


nocontrast        = cfg_const;
nocontrast.tag    = 'nocontrast';
nocontrast.name   = 'No Contrast';
nocontrast.val    = {true};
nocontrast.help   = {'No contrast, just compute fiber counts.'};

contrastsel         = cfg_choice;
contrastsel.tag     = 'contrastsel';
contrastsel.name    = 'Contrast for Tractometry';
contrastsel.values  = {filename4 nocontrast};
contrastsel.help    = {'Select a contrast along which statistics is done.'};



% ---------------------------------------------------------------------
% roinumbers
% ---------------------------------------------------------------------
roinumber        = cfg_entry;
roinumber.tag     = 'roinumber';
roinumber.name    = 'ROI indices';
roinumber.help    = {'Select the ROI indices (mask numbers) used for fiber selection (empty - all ROIs/masks are used).'};
roinumber.strtype = 'e';
roinumber.num     = [0 inf];
roinumber.val = {[]};


% ---------------------------------------------------------------------
% tractStats
% ---------------------------------------------------------------------
tractStat         = cfg_exbranch;
tractStat.tag     = 'tractStats';
tractStat.name    = 'Tract Statistics (Tractometry)';
tractStat.help    = {'help yourself'};
tractStat.val     = {filename1 contrastsel newfilename2};
tractStat.prog    = @(job)dti_ftrop_ui('tractStats',job);
tractStat.vout    = @vout;



% ---------------------------------------------------------------------
% tractStats 
% ---------------------------------------------------------------------
tractStatROI         = cfg_exbranch;
tractStatROI.tag     = 'tractStatsROI';
tractStatROI.name    = 'Tract Statistics with Fiber Selection';
tractStatROI.help    = {'help yourself'};
tractStatROI.val     = {filename1 contrastsel filename5 roinumber seltype fuzzysel newfilename2};
tractStatROI.prog    = @(job)dti_ftrop_ui('tractStatsROI',job);
tractStatROI.vout    = @vout;






% ---------------------------------------------------------------------
% deformFTR 
% ---------------------------------------------------------------------
deformFTR         = cfg_exbranch;
deformFTR.tag     = 'deformFTR';
deformFTR.name    = 'Deformation of Fibers';
deformFTR.help    = {'Non Rigid deformation of fibers.'};
deformFTR.val     = {filename1 filedef filedefi newfilename};
deformFTR.prog    = @(job)dti_ftrop_ui('deformFTR',job);
deformFTR.vout    = @vout;


% ---------------------------------------------------------------------
% mapop 
% ---------------------------------------------------------------------
ftrop         = cfg_choice;
ftrop.tag     = 'ftrop';
ftrop.name    = 'Operations with Streamline Tracts';
ftrop.help    = {'Operations with streamline tracts (Mori Tracts).'};
ftrop.values  = {selectbyROI selectbyOverlapROI eliminatebyROI visitMap visitMapext visitMask endpointMask tractStat tractStatROI deformFTR};

% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
function dep = vout(job)
dep            = cfg_dep;
dep.sname      = 'FTRstruct After Operation';
dep.src_output = substruct('.','files');
dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});




function dep = voutFDext(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'FD (rgb) mrstruct/nifti';
dep(1).src_output = substruct('.','fdrgb');
dep(1).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'FD mrstruct/nifti';
dep(2).src_output = substruct('.','fd');
dep(2).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});

dep(3)            = cfg_dep;
dep(3).sname      = 'EP  mrstruct/nifti';
dep(3).src_output = substruct('.','ep');
dep(3).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});


