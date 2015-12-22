function ex_ftrstruct_Curves2dx = impexp_cfg_ftrstruct_Curves2dx
% 'Export Mori to DX' - MATLABBATCH configuration
% This MATLABBATCH configuration file has been generated automatically
% by MATLABBATCH using ConfGUI. It describes menu structure, validity
% constraints and links to run time code.
% Changes to this file will be overwritten if the ConfGUI batch is executed again.
% Created at 2008-09-02 16:02:44.
% ---------------------------------------------------------------------
% ftrname FTR data file
% ---------------------------------------------------------------------
ftrname         = cfg_files;
ftrname.tag     = 'ftrname';
ftrname.name    = 'FTR data file';
ftrname.filter = 'mat';
ftrname.ufilter = '.*';
ftrname.num     = [1 1];
% ---------------------------------------------------------------------
% ftrstruct FTR data struct
% ---------------------------------------------------------------------
ftrstruct         = cfg_entry;
ftrstruct.tag     = 'ftrstruct';
ftrstruct.name    = 'FTR data struct';
ftrstruct.strtype = 'e';
ftrstruct.num     = [1  1];
% ---------------------------------------------------------------------
% ftr Data
% ---------------------------------------------------------------------
ftr         = cfg_choice;
ftr.tag     = 'ftr';
ftr.name    = 'Data';
ftr.values  = {ftrname ftrstruct };
% ---------------------------------------------------------------------
% index Index
% ---------------------------------------------------------------------
index         = cfg_entry;
index.tag     = 'index';
index.name    = 'Index';
index.strtype = 'n';
index.num     = [1  Inf];
% ---------------------------------------------------------------------
% names Name
% ---------------------------------------------------------------------
names1         = cfg_entry;
names1.tag     = 'names';
names1.name    = 'Name';
names1.strtype = 's';
names1.num     = [1  Inf];
% ---------------------------------------------------------------------
% names Names
% ---------------------------------------------------------------------
names         = cfg_repeat;
names.tag     = 'names';
names.name    = 'Names';
names.values  = {names1 };
names.num     = [1 Inf];
% ---------------------------------------------------------------------
% all All
% ---------------------------------------------------------------------
all         = cfg_const;
all.tag     = 'all';
all.name    = 'All';
all.val     = {1};
% ---------------------------------------------------------------------
% fibers Fiber Selection
% ---------------------------------------------------------------------
fibers         = cfg_choice;
fibers.tag     = 'fibers';
fibers.name    = 'Fiber Selection';
fibers.values  = {index names all };
% ---------------------------------------------------------------------
% ex_ftrstruct_Curves2dx Export Mori to DX
% ---------------------------------------------------------------------
ex_ftrstruct_Curves2dx         = cfg_exbranch;
ex_ftrstruct_Curves2dx.tag     = 'ex_ftrstruct_Curves2dx';
ex_ftrstruct_Curves2dx.name    = 'Export Mori to DX';
ex_ftrstruct_Curves2dx.val     = {ftr fibers };
ex_ftrstruct_Curves2dx.prog = @(job)impexp_run_ftrstruct_Curves2dx(job,'run');
ex_ftrstruct_Curves2dx.vout = @(job)impexp_run_ftrstruct_Curves2dx(job,'vout');
