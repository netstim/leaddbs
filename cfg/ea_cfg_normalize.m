% configure file for calculation of tensors
% DTI Configuration file
% BATCH system (Volkmar Glauche)
%
% File created by Andreas Horn


function normalize = ea_cfg_normalize

% ---------------------------------------------------------------------
% foldername Input raw data filename
% ---------------------------------------------------------------------
foldername         = cfg_files;
foldername.tag     = 'foldername';
foldername.name    = 'Patient folder name';
foldername.help    = {'Select the patient folder name'};
foldername.filter  = 'dir';
foldername.num     = [1 1];


% ---------------------------------------------------------------------
% method Which method is being used for Normalization. 
% ---------------------------------------------------------------------
method         = cfg_menu;
method.tag     = 'method';
method.name    = 'Method';
method.help    = {'Select the method to normalize the nifti images into MNI space'};
method.labels  = {'Schönecker 2009'
                'Schönecker 2009 ? include Preop data'
               'Witt 2013'}';
method.values  = {1 2 3};




% ---------------------------------------------------------------------
% normalize 
% ---------------------------------------------------------------------
normalize         = cfg_exbranch;
normalize.tag     = 'normalize';
normalize.name    = 'Normalize';
normalize.val     = {foldername method};
normalize.help    = {'Normalizes the MR-images to MNI-Space.'};
normalize.prog    = @ea_ui_normalize;
normalize.vout    = @vout;
% ---------------------------------------------------------------------
