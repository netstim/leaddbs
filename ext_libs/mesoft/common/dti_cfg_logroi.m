% configure file for operations on streamline trackings results
%
% for BATCH EDITOR system (Volkmar Glauche)
%
% File created by Susanne Schnell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function logroi = dti_cfg_logroi

% ---------------------------------------------------------------------
% newfilename name of resulting mask after operation
% ---------------------------------------------------------------------
newfilename         = cfg_entry;
newfilename.tag     = 'newfilename';
newfilename.name    = 'New File Name';
newfilename.val     = {'.txt'};
newfilename.help    = {'Type in the name of textfile to be generated.'};
newfilename.strtype = 's';
newfilename.num     = [1 Inf];

% ---------------------------------------------------------------------
% path of log file
% ---------------------------------------------------------------------
logpath      = cfg_files;
logpath.tag  = 'logpath';
logpath.name = 'Output directory';
logpath.help = {'Select the output directory of the textfile.'};
logpath.filter = 'dir';
logpath.num  = [1 1];

% ---------------------------------------------------------------------
% mask, select the correct mask by a number
% ---------------------------------------------------------------------
mask1         = cfg_entry;
mask1.tag     = 'mask1';
mask1.name    = 'mask number';
mask1.help    = {'Select the mask by a number. You need to remember the mask position by yourself.'};
mask1.strtype = 'e';
mask1.num     = [1 Inf];

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
% Log statistics of ROIs
% ---------------------------------------------------------------------
StatsRois         = cfg_branch;
StatsRois.tag     = 'StatsRois';
StatsRois.name    = 'Statistics of all ROIs';
StatsRois.help    = {'Log the statistics of all ROIs.'};



% ---------------------------------------------------------------------
% Log values of one ROI
% ---------------------------------------------------------------------
ValsRoi         = cfg_branch;
ValsRoi.tag     = 'ValsRoi';
ValsRoi.name    = 'Values of ROI';
ValsRoi.help    = {'Log the values of one ROI.'};
ValsRoi.val     = {mask1};


% ---------------------------------------------------------------------
% status
% ---------------------------------------------------------------------
status         = cfg_choice;
status.tag     = 'status';
status.name    = 'Select what to log';
status.help    = {'Select whether you want to log'};
status.values  = {StatsRois  ValsRoi};

% ---------------------------------------------------------------------
% dtdname
% ---------------------------------------------------------------------
dtdname         = cfg_files;
dtdname.tag     = 'dtdname';
dtdname.name    = 'Load DTD';
dtdname.help    = {'Select the corresponding DTD.'};
dtdname.filter  = 'mat';
dtdname.ufilter = '_DTD.*';
dtdname.num     = [1 1];


% ---------------------------------------------------------------------
% Log ROI
% ---------------------------------------------------------------------
logroi        = cfg_exbranch;
logroi.tag     = 'logroi';
logroi.name    = 'Log ROI';
logroi.help    = {'Log the stats or values of all or only one ROI.'};
logroi.val     = {dtdname roiname logpath newfilename status};
logroi.prog    = @(job)dti_logroi_ui(job);
logroi.vout    = @vout;

% ---------------------------------------------------------------------
function dep = vout(job)
dep            = cfg_dep;
dep.sname      = 'MASKstruct After Operation';
dep.src_output = substruct('.','files');
dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});