% configure file for operations on streamline trackings results
%
% for BATCH EDITOR system (Volkmar Glauche)
%
% File created by Susanne Schnell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function unring = dti_cfg_unring


niiname         = cfg_files;
niiname.tag     = 'niiname';
niiname.name    = 'Nifti Files';
niiname.help    = {'Give a series of Niftis/mrstructs to be unringed'};
niiname.filter  = '.nii|.mat';
niiname.ufilter = '.*';
niiname.num     = [1 inf];

% ---------------------------------------------------------------------
% Log ROI
% ---------------------------------------------------------------------
unring        = cfg_exbranch;
unring.tag     = 'realigndef';
unring.name    = 'Unringing';
unring.help    = {'....'};
unring.val     = {niiname};
unring.prog    = @(job)dti_unring(job);
unring.vout    = @vout;

% ---------------------------------------------------------------------

function dep = vout(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'forward deformation';
dep(1).src_output = substruct('.','finames');
dep(1).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});




function out = dti_unring(P)

ringRemoval(P.niiname);
out.finames = P.niiname;
