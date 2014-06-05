% configure file for calculation of tensors
% DTI Configuration file
% BATCH system (Volkmar Glauche)
%
% File created by Andreas Horn


function manualheight = ea_cfg_manualheight

% ---------------------------------------------------------------------
% manualheight Patient folder name
% ---------------------------------------------------------------------
foldername         = cfg_files;
foldername.tag     = 'foldername';
foldername.name    = 'Patient folder';
foldername.help    = {'Select the folder in which the anatomical images are stored.'};
foldername.filter  = 'dir';
foldername.num     = [1 1];


% ---------------------------------------------------------------------
% tensor 
% ---------------------------------------------------------------------
manualheight         = cfg_exbranch;
manualheight.tag     = 'manualheight';
manualheight.name    = 'Manual Height Correction';
manualheight.val     = {foldername};
manualheight.help    = {'Opens a figure to manually change electrode positions.'};
manualheight.prog    = @ea_ui_manualheight;
manualheight.vout    = @vout;