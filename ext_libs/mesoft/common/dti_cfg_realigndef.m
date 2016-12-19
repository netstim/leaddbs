% configure file for operations on streamline trackings results
%
% for BATCH EDITOR system (Volkmar Glauche)
%
% File created by Susanne Schnell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function realigndefi = dti_cfg_realigndef


yname         = cfg_files;
yname.tag     = 'yname';
yname.name    = 'Load forward Deformation';
yname.help    = {'......'};
yname.filter  = 'nii';
yname.ufilter = '.*';
yname.num     = [1 1];

iyname         = cfg_files;
iyname.tag     = 'iyname';
iyname.name    = 'Load inverse Deformation';
iyname.help    = {'......'};
iyname.filter  = 'nii';
iyname.ufilter = '.*';
iyname.num     = [1 1];


matname         = cfg_files;
matname.tag     = 'matname';
matname.name    = 'Matrix for Realignment';
matname.help    = {'......'};
matname.filter  = 'mat';
matname.ufilter = '.*';
matname.num     = [1 1];

% ---------------------------------------------------------------------
% Log ROI
% ---------------------------------------------------------------------
realigndefi        = cfg_exbranch;
realigndefi.tag     = 'realigndef';
realigndefi.name    = 'Realignment of Deformation Field';
realigndefi.help    = {'....'};
realigndefi.val     = {yname iyname matname};
realigndefi.prog    = @(job)dti_realigndef_ui(job);
realigndefi.vout    = @vout;

% ---------------------------------------------------------------------

function dep = vout(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'forward deformation';
dep(1).src_output = substruct('.','finame_y');
dep(1).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'inverse deformation';
dep(2).src_output = substruct('.','finame_iy');
dep(2).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});


