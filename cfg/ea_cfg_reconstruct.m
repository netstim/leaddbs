% configure file for calculation of tensors
% DTI Configuration file
% BATCH system (Volkmar Glauche)
%
% File created by Andreas Horn


function recon = ea_cfg_reconstruct

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
% entrypoint STN/GPi/ViM or Cg25? 
% ---------------------------------------------------------------------
entrypoint         = cfg_menu;
entrypoint.tag     = 'entrypoint';
entrypoint.name    = 'Entrypoint';
entrypoint.val     = {1};
entrypoint.help    = {'Choose which Target region is used to define entry point.'};
entrypoint.labels = {'STN, GPi or ViM'
               'Cg25'}';
entrypoint.values = {1 2};

% ---------------------------------------------------------------------
% Axis contrast 
% ---------------------------------------------------------------------
acontrast         = cfg_menu;
acontrast.tag     = 'acontrast';
acontrast.name    = 'Axis contrast';
acontrast.val     = {4};
acontrast.help    = {'Choose which [combination of] MR-image to use for trajectory reconstruction.'};
acontrast.labels = {'*** Recommended for Axis contrast'
'Use transversal image only'
'Use transversal but smooth'
'Use average of cor & tra, smoothed'
'*** Recommended for z-Contrast'
'Use coronal image only'
'Use coronal but smooth'
'Use average of cor & tra'
'Use cor .* tra'
'Use cor .^4'
'Use (tra .* cor) .^4'
'Use (tra .* cor) .^4 and smooth'
}';
acontrast.values = {1 2 3 4 5 6 7 8 9 10 11 12};

% ---------------------------------------------------------------------
% Coronar contrast 
% ---------------------------------------------------------------------
zcontrast         = cfg_menu;
zcontrast.tag     = 'zcontrast';
zcontrast.name    = 'Z-contrast';
zcontrast.val     = {6};
zcontrast.help    = {'Choose which [combination of] MR-image to use for electrode contact reconstruction.'};
zcontrast.labels = {'*** Recommended for Axis contrast'
'Use transversal image only'
'Use transversal but smooth'
'Use average of cor & tra, smoothed'
'*** Recommended for z-Contrast'
'Use coronal image only'
'Use coronal but smooth'
'Use average of cor & tra'
'Use cor .* tra'
'Use cor .^4'
'Use (tra .* cor) .^4'
'Use (tra .* cor) .^4 and smooth'
}';
zcontrast.values = {1 2 3 4 5 6 7 8 9 10 11 12};

% ---------------------------------------------------------------------
% endtolerance  
% ---------------------------------------------------------------------
endtolerance         = cfg_entry;
endtolerance.tag     = 'endtolerance';
endtolerance.name    = 'Endtolerance';
endtolerance.val     = {10};
endtolerance.help    = {'Enter the number of slices to iterate through, even if information gets fuzzy.','This parameter is seldomly changed, usually, the default of 10 is a good choice.'};
endtolerance.strtype = 'e';
endtolerance.num     = [1 1];

% ---------------------------------------------------------------------
% Distance Threshold  
% ---------------------------------------------------------------------
distancethresh         = cfg_entry;
distancethresh.tag     = 'distancethresh';
distancethresh.name    = 'Distance Threshold';
distancethresh.val     = {4};
distancethresh.help    = {'A parameter affecting the stop criterion.','This parameter is seldomly changed, usually, the default of 10 is a good choice.'};
distancethresh.strtype = 'e';
distancethresh.num     = [1 1];

% ---------------------------------------------------------------------
% Transv. Std. Cutoff  
% ---------------------------------------------------------------------
tstdcut         = cfg_entry;
tstdcut.tag     = 'tstdcut';
tstdcut.name    = 'Transv. Std-cutoff.';
tstdcut.val     = {0.9};
tstdcut.help    = {'A parameter affecting the thresholding of each slice.','This parameter is seldomly changed, usually, the default of 10 is a good choice.'};
tstdcut.strtype = 'e';
tstdcut.num     = [1 1];

% ---------------------------------------------------------------------
% Coron. Std. Cutoff  
% ---------------------------------------------------------------------
cstdcut         = cfg_entry;
cstdcut.tag     = 'cstdcut';
cstdcut.name    = 'Cor. Std-cutoff.';
cstdcut.val     = {1.0};
cstdcut.help    = {'A parameter affecting the height reconstruction of the electrodes.','This parameter is seldomly changed, usually, the default of 10 is a good choice.'};
cstdcut.strtype = 'e';
cstdcut.num     = [1 1];

% ---------------------------------------------------------------------
% recon 
% ---------------------------------------------------------------------
recon         = cfg_exbranch;
recon.tag     = 'recon';
recon.name    = 'Do Reconstruction';
recon.val     = {foldername entrypoint acontrast zcontrast endtolerance distancethresh tstdcut cstdcut};
recon.help    = {'Calculates the DTI tensors.'};
recon.prog    = @ea_ui_reconstruct;
recon.vout    = @vout;
% ---------------------------------------------------------------------

