
% configure file for calculation of tensors
% DTI Configuration file
% BATCH system (Volkmar Glauche)
%
% File created by Susanne Schnell


function tracking = dti_cfg_tracking

% ---------------------------------------------------------------------
% filename Input dtdstruct
% ---------------------------------------------------------------------
filename         = cfg_files;
filename.tag     = 'filename';
filename.name    = 'DTD file name';
filename.help    = {'Select the calculated tensor file from the file selection menue("_DTD.mat").'};
filename.filter  = 'mat';
filename.ufilter = '_DTD\.mat$';
filename.num     = [1 1];

% ---------------------------------------------------------------------
% filename Input FTRstruct
% ---------------------------------------------------------------------
filenameFTR         = cfg_files;
filenameFTR.tag     = 'filenameFTR';
filenameFTR.name    = 'FTR file name';
filenameFTR.help    = {'Select the GT-FTR.'};
filenameFTR.filter  = 'mat';
filenameFTR.ufilter = '_FTR\.mat$';
filenameFTR.num     = [1 1];


% ---------------------------------------------------------------------
% startname, select the correct mask by name
% ---------------------------------------------------------------------
startname         = cfg_entry;
startname.tag     = 'startname';
startname.name    = 'Type in the name of the mask';
startname.help    = {'Type in the name of the mask.'};
startname.strtype = 's';
startname.num     = [1 Inf];


% ---------------------------------------------------------------------
% startmask, select the correct mask by a number
% ---------------------------------------------------------------------
startnumber         = cfg_entry;
startnumber.tag     = 'startnumber';
startnumber.name    = 'Mask Number';
startnumber.help    = {['Select the mask by a number.'...
                    'You need to remember the mask position by yourself.']};
startnumber.strtype = 'e';
startnumber.num     = [1 1];

% ---------------------------------------------------------------------
% startmask, select the correct mask by a number
% ---------------------------------------------------------------------
startmask         = cfg_choice;
startmask.tag     = 'startmask';
startmask.name    = 'Select the mask by number or by name';
startmask.help    = {['Select the mask either by a number or by the name.' ...
                        'The name has to be spelled correctly. For the number'...
                        'you need to remember the mask position by yourself.']};
startmask.values = {startname startnumber};

% ---------------------------------------------------------------------
% startfileROI as start mask (mori)
% ---------------------------------------------------------------------
startfile         = cfg_files;
startfile.tag     = 'startfile';
startfile.name    = 'Load ROI as Start Mask';
startfile.help    = {'Select a maskstruct as start mask, can be some region of interest.'};
startfile.filter  = 'mat';
startfile.ufilter = '.*';
startfile.num     = [1 1];

% ---------------------------------------------------------------------
% startdef starting criteria as mask (mori)
% ---------------------------------------------------------------------
startdef         = cfg_branch;
startdef.tag     = 'startdef';
startdef.name    = 'Definition of Start Region';
startdef.help    = {'Chose the mask-file as starting criteria and select the mask.'};
startdef.val     = {startfile startmask};

% ---------------------------------------------------------------------
% stopname, select the correct mask by name
% ---------------------------------------------------------------------
stopname         = cfg_entry;
stopname.tag     = 'stopname';
stopname.name    = 'Type in the name of the mask';
stopname.help    = {'Type in the name of the mask.'};
stopname.strtype = 's';
stopname.num     = [1 Inf];

% ---------------------------------------------------------------------
% stopnumber, select the correct mask by a number
% ---------------------------------------------------------------------
stopnumber         = cfg_entry;
stopnumber.tag     = 'stopnumber';
stopnumber.name    = 'Mask Number';
stopnumber.help    = {'Select the mask by a number.'};
stopnumber.strtype = 'e';
stopnumber.num     = [1 1];

% ---------------------------------------------------------------------
% stopmask, select the correct mask by a number
% ---------------------------------------------------------------------
stopmask         = cfg_choice;
stopmask.tag     = 'stopmask';
stopmask.name    = 'Mask Number or Name';
stopmask.help    = {'Select the mask by a number or by name.'};
stopmask.values  = {stopnumber stopname};

% ---------------------------------------------------------------------
% stopmask ROI as stop mask (mori)
% ---------------------------------------------------------------------
stopfile         = cfg_files;
stopfile.tag     = 'stopfile';
stopfile.name    = 'Load ROI as Stop Mask';
stopfile.help    = {'Load a file containing the stop mask.'};
stopfile.filter  = 'mat';
stopfile.ufilter = '.*';
stopfile.num     = [1 1];

% ---------------------------------------------------------------------
% stopdef stopping criteria as mask (mori)
% ---------------------------------------------------------------------
stopdef         = cfg_branch;
stopdef.tag     = 'stopdef';
stopdef.name    = 'Definition of Stop Region';
stopdef.help    = {'Chose the mask-file as stopping criteria and select the mask.'};
stopdef.val     = {stopfile stopmask};

% ---------------------------------------------------------------------
% seedname, select the correct mask by name
% ---------------------------------------------------------------------
seedname         = cfg_entry;
seedname.tag     = 'seedname';
seedname.name    = 'Type in the name of the mask';
seedname.help    = {'Type in the name of the mask from which you want to start tracking.'};
seedname.strtype = 's';
seedname.num     = [1 Inf];

% ---------------------------------------------------------------------
% seedmask, select the correct mask by a number (prob)
% ---------------------------------------------------------------------
seednumber         = cfg_entry;
seednumber.tag     = 'seednumber';
seednumber.name    = 'Mask Number';
seednumber.help    = {'Select the mask by a number. You need to remember the position number by yourself!', ...
                      'If more then one number is given, tracking will start from each of the specified masks.',...
                      'Specify ''Inf'' to use all masks in the given mask struct.'};
seednumber.strtype = 'e';
seednumber.num     = [1 Inf];

% ---------------------------------------------------------------------
% seedmask, select the correct mask by a number (prob)
% ---------------------------------------------------------------------
seedmask         = cfg_choice;
seedmask.tag     = 'seedmask';
seedmask.name    = 'Mask Number or Name';
seedmask.help    = {'Select the mask by a number or by name. You need to remember the correct spelling. In case of'...
                    'a number you need to remember the position number by yourself!', ...
                    'If more then one number is given, tracking will start from each of the specified masks.',...
                    'Specify ''Inf'' to use all masks in the given mask struct.'};
seedmask.values  = {seednumber seedname};

% ---------------------------------------------------------------------
% seedfile, filename of seedroi (prob)
% ---------------------------------------------------------------------
seedfile         = cfg_files;
seedfile.tag     = 'seedfile';
seedfile.name    = 'Mask filename';
seedfile.help    = {'Select a mask specifying the seed region.'};
seedfile.filter  = 'mat';
seedfile.ufilter = '.*';
seedfile.num     = [1 1];

% ---------------------------------------------------------------------
% seedroi MaskStruct containing Start Region (prob)
% ---------------------------------------------------------------------
seedroi         = cfg_branch;
seedroi.tag     = 'seedroi';
seedroi.name    = 'MaskStruct containing Start Region';
seedroi.help    = {'Chose the mask-file for the definition of the seed region and select the mask.'};
seedroi.val     = {seedfile seedmask};

% ---------------------------------------------------------------------
% seedposition enter position of seed voxel (prob)
% ---------------------------------------------------------------------
seedposition         = cfg_entry;
seedposition.tag     = 'seedposition';
seedposition.name    = 'Enter Position of the Seed Voxel';
seedposition.val{1}  = [int16(1) int16(1) int16(1)];
seedposition.help    = {'Chose the mask-file for the definition of the seed region and select the mask.'};
seedposition.strtype = 'e';
seedposition.num     = [1 3];

% ---------------------------------------------------------------------
% seedchoice Definition of Start Region (prob)
% ---------------------------------------------------------------------
seedchoice      = cfg_choice;
seedchoice.tag  = 'seedchoice';
seedchoice.name = 'Definition of Start Region';
seedchoice.help = {['The seed region can be specified in different ' ...
                    'ways.']};
seedchoice.values = {seedroi seedposition};

% ---------------------------------------------------------------------
% defname, select the correct mask by name
% ---------------------------------------------------------------------
defname         = cfg_entry;
defname.tag     = 'defname';
defname.name    = 'Type in the name of the mask';
defname.help    = {'Type in the name of the mask in which tracking should be performed.'};
defname.strtype = 's';
defname.num     = [1 Inf];


% ---------------------------------------------------------------------
% defmask, select the correct mask by a number (prob)
% ---------------------------------------------------------------------
defnumber         = cfg_entry;
defnumber.tag     = 'defnumber';
defnumber.name    = 'Mask Number';
defnumber.help    = {'Select the mask by a number. You need to remember the correct spelling. In case of',...
                     'a number you need to remember the position number by yourself!',...
                     'Or if you created and saved the mask in the fiberviewer_tool, check at which position',...
                     'or under which name your mask is listed.'};
defnumber.strtype = 'e';
defnumber.num     = [1 1];

% ---------------------------------------------------------------------
% defmask, select the correct mask by a number (prob)
% ---------------------------------------------------------------------
defmask         = cfg_choice;
defmask.tag     = 'defmask';
defmask.name    = 'Mask Number or Name';
defmask.help    = {'Select the mask by a number or by name.  You need to remember the position number by yourself! If you created and saved the mask in the fiberviewer_tool, check at which position your mask is listed.'};
defmask.values  = {defnumber defname};


% ---------------------------------------------------------------------
% deffile definition of tracking area (saved in mask) (prob)
% ---------------------------------------------------------------------
deffile         = cfg_files;
deffile.tag     = 'deffile';
deffile.name    = 'Mask filename';
deffile.help    = {'Select the tracking mask in the file selection dialogue.'};
deffile.filter  = 'mat';
deffile.ufilter = '.*';
deffile.num     = [1 1];

% ---------------------------------------------------------------------
% defarea Input tracking area (prob)
% ---------------------------------------------------------------------
defarea         = cfg_branch;
defarea.tag     = 'defarea';
defarea.name    = 'Definition of Tracking Area using a Mask';
defarea.help    = {'Chose a mask file and select the mask defining the tracking area (e.g. white matter mask)'};
defarea.val     = {deffile defmask};


% ---------------------------------------------------------------------
% deffa, definition of tracking area using FA (prob)
% ---------------------------------------------------------------------
deffa         = cfg_entry;
deffa.tag     = 'deffa';
deffa.name    = 'FA threshold';
deffa.val     = {0.1};
deffa.help    = {'Enter the FA value above which tracking will be performed only. default = 0.1'};
deffa.strtype = 'e';
deffa.num     = [1 1];

% ---------------------------------------------------------------------
% deffa, definition of tracking area using Trace (prob)
% ---------------------------------------------------------------------
deftrace         = cfg_entry;
deftrace.tag     = 'deftrace';
deftrace.name    = 'Trace threshold';
deftrace.val     = {0.002};
deftrace.help    = {'Enter the Trace (mean diffusivity) value below which tracking will be performed only. default = 0.002'};
deftrace.strtype = 'e';
deftrace.num     = [1 1];

% ---------------------------------------------------------------------
% defarea Input tracking area (prob)
% ---------------------------------------------------------------------
defthresh         = cfg_branch;
defthresh.tag     = 'defthresh';
defthresh.name    = 'Definition of Tracking Area using Thresholds';
defthresh.help    = {'Select thresholds for FA (tracking only above this threshold) and Trace (tracking only below this threshold).'};
defthresh.val     = {deffa deftrace};

% ---------------------------------------------------------------------
% start Which criteria for starting tracking, mask or thresholds? (prob)
% ---------------------------------------------------------------------
defroi         = cfg_choice;
defroi.name    = 'Define Tracking Area';
defroi.tag     = 'defroi';
defroi.values  = {defarea defthresh};
defroi.help    = {'Do you want to define the tracking area by thresholds or by masks?'};


% ---------------------------------------------------------------------
% fname name of resulting file containing the fibre tracks (prob)
% ---------------------------------------------------------------------
fname         = cfg_entry;
fname.tag     = 'fname';
fname.name    = 'File Name of the resulting tracks';
fname.help    = {'Type in the name of the new file containg the fibre tracks.'};
fname.strtype = 's';
fname.num     = [1 Inf];

% ---------------------------------------------------------------------
% dir Output directory (prob)
% ---------------------------------------------------------------------
dir      = cfg_files;
dir.tag  = 'dir';
dir.name = 'Output directory';
dir.help = {'Select the output directory.'};
dir.filter = 'dir';
dir.num  = [1 1];

% ---------------------------------------------------------------------
% out User-specified output location (prob)
% ---------------------------------------------------------------------
out      = cfg_branch;
out.tag  = 'out';
out.name = 'User-specified output location';
out.val  = {dir fname};
out.help = {'Specify output directory and filename.'};

% ---------------------------------------------------------------------
% auto Automatically determine output filename (prob)
% ---------------------------------------------------------------------
auto      = cfg_const;
auto.tag  = 'auto';
auto.name = 'Automatically determine output filename';
auto.val  = {1};

% ---------------------------------------------------------------------
% newprobfile Select output location (prob)
% ---------------------------------------------------------------------
newmorifile        = cfg_choice;
newmorifile.tag    = 'newmorifile';
newmorifile.name   = 'Select output location';
newmorifile.values = {auto out};
newmorifile.help   = {'Name and location can be specified explicitly, or they can be derived from the input DTD struct filename.'};

% ---------------------------------------------------------------------
% newprobfile Select output location (prob)
% ---------------------------------------------------------------------
newprobfile        = cfg_choice;
newprobfile.tag    = 'newprobfile';
newprobfile.name   = 'Select output location';
newprobfile.values = {auto out};
newprobfile.help   = {'Name and location can be specified explicitly, or they can be derived from the input DTD struct filename.'};

% ---------------------------------------------------------------------
% maxangle maximum allowed angle of fibre track (mori)
% ---------------------------------------------------------------------
maxangle         = cfg_entry;
maxangle.tag     = 'maxangle';
maxangle.name    = 'Maximum Allowed Angle';
maxangle.val     = {53.1};
maxangle.help    = {'Type in the maximum angle a fibre tract can bend. default = 53.1'};
maxangle.strtype = 'e';
maxangle.num     = [1 1];

% ---------------------------------------------------------------------
% minvox minumum amount of voxels for fibre tracts (mori)
% ---------------------------------------------------------------------
minvox         = cfg_entry;
minvox.tag     = 'minvox';
minvox.name    = 'Minimum amount of voxels';
minvox.val     = {5};
minvox.help    = {'Enter the minimum length of the tract in amount of voxels. default = 5'};
minvox.strtype = 'e';
minvox.num     = [1 1];

% ---------------------------------------------------------------------
% startfaLim start tracking if voxel has an FA value below (mori)
% ---------------------------------------------------------------------
startfaLim         = cfg_entry;
startfaLim.tag     = 'startfaLim';
startfaLim.name    = 'Start at FA Limit';
startfaLim.val     = {0.25};
startfaLim.help    = {'Stop if voxel has an FA value above. default = 0.25'};
startfaLim.strtype = 'e';
startfaLim.num     = [1 1];

% ---------------------------------------------------------------------
% startTrLim start tracking if voxel has a Trace value below (mori)
% ---------------------------------------------------------------------
startTrLim         = cfg_entry;
startTrLim.tag     = 'startTrLim';
startTrLim.name    = 'Start at Trace Limit';
startTrLim.val     = {0.0016};
startTrLim.help    = {'Start if voxel has an Trace (= mean diffusivity) value below. default = 0.0016'};
startTrLim.strtype = 'e';
startTrLim.num     = [1 1];


% ---------------------------------------------------------------------
% startthreshold select start criteria via thresholds (mori)
% ---------------------------------------------------------------------
startthreshold         = cfg_branch;
startthreshold.tag     = 'startthreshold';
startthreshold.name    = 'Start Criteria via Threshold';
startthreshold.help    = {'Select FA and Trace(D) thresholds as starting criteria for the streamline tracking.'};
startthreshold.val     = {startfaLim startTrLim};

% ---------------------------------------------------------------------
% start Which criteria for starting tracking, mask or thresholds? (mori)
% ---------------------------------------------------------------------
start         = cfg_choice;
start.name    = 'Select Start Criteria';
start.tag     = 'start';
start.values  = {startdef startthreshold};
start.help    = {'Chose between a mask or FA and Trace(D) thresholds as starting criteria.'};

% ---------------------------------------------------------------------
% stopfaLim stop tracking if voxel has an FA value below (mori)
% ---------------------------------------------------------------------
stopfaLim         = cfg_entry;
stopfaLim.tag     = 'stopfaLim';
stopfaLim.name    = 'Stop at FA Limit';
stopfaLim.val     = {0.15};
stopfaLim.help    = {'Stop if voxel has an FA value below. default = 0.15'};
stopfaLim.strtype = 'e';
stopfaLim.num     = [1 1];

% ---------------------------------------------------------------------
% stopTrLim stop tracking if voxel has a Trace value below (mori)
% ---------------------------------------------------------------------
stopTrLim         = cfg_entry;
stopTrLim.tag     = 'stopTrLim';
stopTrLim.name    = 'Stop at Trace Limit';
stopTrLim.val     = {0.002};
stopTrLim.help    = {'Stop if voxel has an Trace (= mean diffusivity) value below. default = 0.002'};
stopTrLim.strtype = 'e';
stopTrLim.num     = [1 1];

% ---------------------------------------------------------------------
% startthreshold select start criteria via thresholds (mori)
% ---------------------------------------------------------------------
stopthreshold         = cfg_branch;
stopthreshold.tag     = 'stopthreshold';
stopthreshold.name    = 'Stop Criteria via Thresholds';
stopthreshold.help    = {'Select FA and Trace(D) thresholds as stopping criteria for the streamline tracking.'};
stopthreshold.val     = {stopfaLim stopTrLim};

% ---------------------------------------------------------------------
% stop Which criteria for stopping tracking, mask or thresholds? (mori)
% ---------------------------------------------------------------------
stop         = cfg_choice;
stop.name    = 'Select Stop Criteria';
stop.tag     = 'stop';
stop.values  = {stopdef stopthreshold};
stop.help    = {'Chose between a mask or FA and Trace thresholds as start criteria.'};

% ---------------------------------------------------------------------
% randsampno how many streamlines to start (only if randsampler on yes)
% (mori)
% ---------------------------------------------------------------------
randsampno         = cfg_entry;
randsampno.tag     = 'randsampno';
randsampno.name    = 'How many streamlines per voxel?';
randsampno.val     = {2};
randsampno.help    = {'How many streamlines should be started per voxel? This option is only possible if the Random Sampler is switched on. default = 2'};
randsampno.strtype = 'e';
randsampno.num     = [1 1];


% ---------------------------------------------------------------------
% randsampler random start position within voxel? (mori)
% ---------------------------------------------------------------------
randsampler         = cfg_menu;
randsampler.tag     = 'randsampler';
randsampler.name    = 'Switch on Random Sampler?';
randsampler.val     = {2};
randsampler.help    = {'Do you want more than one streamline, which start in random position in the voxel or not? default: No'};
randsampler.labels  = {'Yes - start at random position in voxels'
               'No -  start in the middle of voxels'}';
randsampler.values  = {1 2};


% ---------------------------------------------------------------------
% threshold threshold in B0 images at which tensors are defined as 0 (prob)
% ---------------------------------------------------------------------
exponent         = cfg_entry;
exponent.tag     = 'exponent';
exponent.name    = 'Exponent';
exponent.val     = {4};
exponent.help    = {'Enter the exponent by which the underlying diffusion tensor should be sharpened. default = 4'};
exponent.strtype = 'e';
exponent.num     = [1 1];

% ---------------------------------------------------------------------
% nowalks number of walks (prob)
% ---------------------------------------------------------------------
nowalks         = cfg_entry;
nowalks.tag     = 'nowalks';
nowalks.name    = 'Number of Walks';
nowalks.val     = {1000};
nowalks.help    = {'Enter the amount of random walks. The default value of 1000 is rather a quick and dirt value. default = 1000'};
nowalks.strtype = 'e';
nowalks.num     = [1 1];

% ---------------------------------------------------------------------
% fibrelength definitin of fibrelength (prob)
% ---------------------------------------------------------------------
fibrelength         = cfg_entry;
fibrelength.tag     = 'fibrelength';
fibrelength.name    = 'Maximum Number of Voxels per Fibre';
fibrelength.val     = {150};
fibrelength.help    = {'Enter the maximum number of voxels a fibre can have. default = 150'};
fibrelength.strtype = 'e';
fibrelength.num     = [1 1];

% ---------------------------------------------------------------------
% probrand chose between two probabilistic algorithms (prob)
% ---------------------------------------------------------------------
probrand         = cfg_menu;
probrand.tag     = 'probrand';
probrand.name    = 'Probabilistic Algorithm';
probrand.val     = {2};
probrand.help    = {'Chose between two probabilistic algorithms both based on the diffusion tensor.'
    'default: 2'};
probrand.labels = {'Probabilistic tracking'
               'Extended probabilistic tracking'}';
probrand.values = {1 2};

% ---------------------------------------------------------------------
% revisits allow looping revisits of a voxel or not (prob)
% ---------------------------------------------------------------------
revisits         = cfg_menu;
revisits.tag     = 'revisits';
revisits.name    = 'Allow Revisits';
revisits.val     = {1};
revisits.help    = {'Select if you want to allow looping revisits of a voxel or not. default: Yes'};
revisits.labels  = {'Yes - allow revisits'
               'No -  do not allow revisits'}';
revisits.values  = {1 0};

% ---------------------------------------------------------------------
% probabilistic Selection of parameters (prob)
% ---------------------------------------------------------------------
probabilistic         = cfg_exbranch;
probabilistic.tag     = 'probabilistic';
probabilistic.name    = 'Probabilistic Tracking';
probabilistic.help    = {'Select the tensor file and additional information necessary for the algorithm.'};
probabilistic.val     = {filename seedchoice defroi nowalks revisits fibrelength probrand exponent newprobfile};
probabilistic.prog    = @dti_tracking_probabilistic_ui;
probabilistic.vout    = @vout;

% ---------------------------------------------------------------------
% mori Input mrstruct file name (mori)
% ---------------------------------------------------------------------
mori         = cfg_exbranch;
mori.tag     = 'mori';
mori.name    = 'Mori Streamline Tracking';
mori.help    = {'Select the dtdstruct file and the enter the tracking related parameters.'};
mori.val     = {filename start stop maxangle randsampler minvox randsampno newmorifile};
mori.prog    = @dti_tracking_mori_ui;
mori.vout    = @vout;

% ---------------------------------------------------------------------
% GT Input mrstruct file name (global tracker)
% ---------------------------------------------------------------------
gtrack         = cfg_exbranch;
gtrack.tag     = 'GTtrack';
gtrack.name    = 'Global Tracking';
gtrack.help    = {'Select the dtdstruct file and the enter the tracking related parameters.'};
gtrack.val     = GTvalList;
gtrack.prog    = @dti_tracking_global_ui;
gtrack.vout    = @voutGT;

% ---------------------------------------------------------------------
% GT fiber accum
% ---------------------------------------------------------------------



numsamps         = cfg_entry;
numsamps.tag     = 'numsamps';
numsamps.name    = 'Number of Samples';
numsamps.val     = {10};
numsamps.help    = {'Enter the number of samples.'};
numsamps.strtype = 'e';
numsamps.num     = [1 1];


numits         = cfg_entry;
numits.tag     = 'numits';
numits.name    = 'Number of Iterations';
numits.val     = {10^6};
numits.help    = {'Enter the number of iterations per sample.'};
numits.strtype = 'e';
numits.num     = [1 1];

temp         = cfg_entry;
temp.tag     = 'temp';
temp.name    = 'Temperature';
temp.val     = {0.1};
temp.help    = {'Temperature where we are sampling.'};
temp.strtype = 'e';
temp.num     = [1 1];

gtrack_accum         = cfg_exbranch;
gtrack_accum.tag     = 'GTtrack_accum';
gtrack_accum.name    = 'Global Tracking / Accumulate FTRs ';
gtrack_accum.help    = {'Select the GT FTR and enter the tracking related parameters.'};
gtrack_accum.val     = {filenameFTR fname numsamps numits temp};
gtrack_accum.prog    = @dti_tracking_global_ui;
takefirst = @(x) x(1);
gtrack_accum.vout    = @(x) takefirst(voutGT(x));




% ---------------------------------------------------------------------
% tracking 
% ---------------------------------------------------------------------
tracking         = cfg_choice;
tracking.tag     = 'tracking';
tracking.name    = 'Fibre Tracking';
tracking.values  = {probabilistic mori gtrack gtrack_accum};
tracking.help    = {'Select the fibre tracking algorithm, select parameters and start tracking.'};
% -------------------------


% ---------------------------------------------------------------------
function dep = vout(job)
% Common to both algorithms
dep            = cfg_dep;
dep.sname      = 'Fibre Tracking';
dep.src_output = substruct('.','files');
dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});

function dep = voutGT(job)
% Common to both algorithms
dep(1)            = cfg_dep;
dep(1).sname      = 'Fibre Tracking (FTR)';
dep(1).src_output = substruct('.','ftr');
dep(1).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'Fibre Tracking (FD)';
dep(2).src_output = substruct('.','fd');
dep(2).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});

dep(3)            = cfg_dep;
dep(3).sname      = 'Fibre Tracking (EPD)';
dep(3).src_output = substruct('.','epd');
dep(3).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});




%%%%%%%%%%%%% GLOBAL TRACKING batch,  mrc 2010


function list = GTvalList

%%%%%%%%%%%%% DTD/HARDI file selection

fnameDTD         = cfg_files;
fnameDTD.tag     = 'filenameDTD';
fnameDTD.name    = 'DTD filename';
fnameDTD.help    = {'Select the calculated tensor file from the file selection menue("_DTD.mat").'};
fnameDTD.filter  = 'mat';
fnameDTD.ufilter = '_DTD\.mat$';
fnameDTD.num     = [1 1];

fnameHARDI         = cfg_files;
fnameHARDI.tag     = 'filenameHARDI';
fnameHARDI.name    = 'HARDI filename';
fnameHARDI.help    = {'Select the HARDI file from the file selection menue("_HARDI.mat").'};
fnameHARDI.filter  = 'mat';
fnameHARDI.ufilter = '_HARDI\.mat$';
fnameHARDI.num     = [1 1];

fnameNII         = cfg_files;
fnameNII.tag     = 'filenameNII';
fnameNII.name    = 'Nifti DW filename';
fnameNII.help    = {'Select Diffusion-weighted Nifti images.'};
fnameNII.filter  = 'image';
fnameNII.ufilter = '.*';
fnameNII.num     = [1 inf];


fname         = cfg_choice;
fname.tag     = 'fname';
fname.name    = 'HARDI/Nifti/DTD filename';
fname.help    = {'Select HARDI,DTD of Nifti for tracking.'};
fname.values     = {fnameDTD fnameNII fnameHARDI};


fnameftr         = cfg_entry;
fnameftr.tag     = 'fname';
fnameftr.name    = 'File Name of the resulting tracks';
fnameftr.help    = {'Type in the name of the new file containg the fibre tracks.'};
fnameftr.strtype = 's';
fnameftr.num     = [1 Inf];

dir      = cfg_files;
dir.tag  = 'dir';
dir.name = 'Output directory';
dir.help = {'Select the output directory.'};
dir.filter = 'dir';
dir.num  = [1 1];

out      = cfg_branch;
out.tag  = 'out';
out.name = 'User-specified output location';
out.val  = {dir fnameftr};
out.help = {'Specify output directory and filename.'};

auto      = cfg_const;
auto.tag  = 'auto';
auto.name = 'Automatically determine output filename';
auto.val  = {1};

newfile        = cfg_choice;
newfile.tag    = 'newfile';
newfile.name   = 'Select output location';
newfile.values = {auto out};
newfile.help   = {'Name and location can be specified explicitly, or they can be derived from the input DTD struct filename.'};



%%%%%%%%%%%%% parameters

parameters         = cfg_menu;
parameters.tag     = 'parameters';
parameters.name    = 'Parameter Suggestion';
parameters.val     = {0};
parameters.help    = {'Decide for a sparse (fast) or dense (slow and accurate) setting'};
parameters.labels = {'Sparse' 'Dense'}';
parameters.values = {0 1};


roiid         = cfg_entry;
roiid.tag     = 'roiid';
roiid.name    = 'name';
roiid.val     = {1};
roiid.help    = {'Each ROI inside a maskstruct has unique name. Enter here the name of the mask which should be used as tracking area.'};
roiid.strtype = 's';
roiid.num     = [1 inf];

%%%%%%%%%%%%% mask selection


estimate         = cfg_const;
estimate.name = 'Estimate Mask';
estimate.tag = 'estimate';
estimate.help = {'A whole brain mask is automatically estimated.'};
estimate.val = {0};

fnameMASK         = cfg_files;
fnameMASK.tag     = 'filenameMASK';
fnameMASK.name    = 'Maskstruct filename';
fnameMASK.help    = {'Select a predetermined tracking area from a maskstruct file (".mat").'};
fnameMASK.filter  = 'mat';
fnameMASK.ufilter = '.mat$';
fnameMASK.num     = [1 1];

fnameMASKnii         = cfg_files;
fnameMASKnii.tag     = 'filenameMASKnii';
fnameMASKnii.name    = 'Nifti filename';
fnameMASKnii.help    = {'Select a predetermined tracking area from a Nifti file.'};
fnameMASKnii.filter  = 'image';
fnameMASKnii.ufilter = '.*';
fnameMASKnii.num     = [1 1];



roiid         = cfg_entry;
roiid.tag     = 'roiid';
roiid.name    = 'name';
roiid.val     = {'untitled'};
roiid.help    = {'Each ROI inside a maskstruct has a unique name. Enter here the name of the mask which should be used as tracking area.'};
roiid.strtype = 's';
roiid.num     = [1 inf];

thresholdMask         = cfg_entry;
thresholdMask.tag     = 'thresholdMask';
thresholdMask.name    = 'Select Threshold';
thresholdMask.val     = {0.5};
thresholdMask.help    = {'Select a threshold (if mask is not binary).'};
thresholdMask.strtype = 'e';
thresholdMask.num     = [1 1];



mask         = cfg_branch;
mask.tag     = 'maskstruct';
mask.name    = 'Tracking Area from maskstruct-file';
mask.help    = {'Select Tracking Area from maskstruct-file.'};
mask.val     = {fnameMASK roiid};

masknii         = cfg_branch;
masknii.tag     = 'masknii';
masknii.name    = 'Tracking Area from Nifti-file';
masknii.help    = {'Select Tracking Area from Nifti-file.'};
masknii.val     = {fnameMASKnii thresholdMask};

trackingarea         = cfg_choice;
trackingarea.tag     = 'trackingarea';
trackingarea.name    = 'Tracking Area';
trackingarea.help    = {'Select Tracking Area.'};
trackingarea.values     = {mask masknii estimate};

%%%%%%%%% custom parameters


auto_para_weight     = cfg_const;
auto_para_weight.name = 'Automatically determined weight based on suggestion';
auto_para_weight.tag = 'auto_para_weight';
auto_para_weight.help = {'The weight parameter is automatically determined during parameter suggestion (not suitable for group studies).'};
auto_para_weight.val = {0};

custom_para_weight   = cfg_entry;
custom_para_weight.tag     = 'custom_para_weight';
custom_para_weight.name    = 'Custom Parameter';
custom_para_weight.val     = {0.05};
custom_para_weight.help    = {'Enter the segment weight'};
custom_para_weight.strtype = 'e';
custom_para_weight.num     = [1 1];

para_weight         = cfg_choice;
para_weight.tag     = 'para_weight';
para_weight.name    = 'Segment Weight';
para_weight.values   = { auto_para_weight custom_para_weight};
para_weight.help    = {'Choose automatically determined weight (not suitable for group studies) or customize the weight parameter'};




auto_para_other     = cfg_const;
auto_para_other.name = 'Automatically determined parameters based on suggestion';
auto_para_other.tag = 'auto_para_other';
auto_para_other.help = {'All other parameters are automatically determined during parameter suggestion.'};
auto_para_other.val = {0};

custom_para_other   = cfg_entry;
custom_para_other.tag     = 'custom_para_other';
custom_para_other.name    = 'Custom Parameters';
custom_para_other.val     = {[0.1, 0.001 , 50 , 5*10^8 , 1 , 2 , 0.2 , 1, 0.5]};
custom_para_other.help    = {'Enter a list of parameteres [start Temp , stop Temp, #Steps, #Iterations, Width(mm), Length(mm), Dens. Penalty, b-value'};
custom_para_other.strtype = 'e';
custom_para_other.num     = [9 1];

para_other         = cfg_choice;
para_other.tag     = 'para_other';
para_other.name    = 'Other Parameters';
para_other.values   = { auto_para_other custom_para_other};
para_other.help    = {'Choose automatically determined parameters or customize them'};





minlen         = cfg_entry;
minlen.tag     = 'minlen';
minlen.name    = 'Minimum fiber length';
minlen.val     = {10};
minlen.help    = {'Enter the minimum number of segments per fiber (a number between 1 and 256)'};
minlen.strtype = 'e';
minlen.num     = [1 1];


maxlen         = cfg_entry;
maxlen.tag     = 'maxlen';
maxlen.name    = 'Maximum fiber length';
maxlen.val     = {inf};
maxlen.help    = {'Enter the maximum number of segments per fiber (a number between 1 and Inf)'};
maxlen.strtype = 'e';
maxlen.num     = [1 1];





list = {fname newfile trackingarea parameters para_weight para_other minlen maxlen};
