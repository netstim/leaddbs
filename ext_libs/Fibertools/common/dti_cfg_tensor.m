% configure file for calculation of tensors
% DTI Configuration file
% BATCH system (Volkmar Glauche)
%
% File created by Susanne Schnell


function tensor = dti_cfg_tensor

% ---------------------------------------------------------------------
% filename Input raw data filename
% ---------------------------------------------------------------------
filename         = cfg_files;
filename.tag     = 'filename';
filename.name    = 'Raw data file name';
filename.help    = {'These are the raw images that will be needed for tensor calculation.'
    'Here only data of type "_raw.bin" or matlab structure "mrstruct" are possible, meaning before tensor calculation you need to run "Read Data"'};
filename.filter  = 'any';
filename.ufilter = '_raw\.bin$';
filename.num     = [1 1];


% ---------------------------------------------------------------------
% threshold threshold in B0 images at which tensors are defined as 0 
% ---------------------------------------------------------------------
threshold         = cfg_entry;
threshold.tag     = 'threshold';
threshold.name    = 'Threshold';
threshold.val     = {40};
threshold.help    = {'Enter the threshold in B0 images at which tensors are defined as 0.'...
                    'default = 40'};
threshold.strtype = 'e';
threshold.num     = [1 1];

% ---------------------------------------------------------------------
% singleslice Save single slice as matlab structure? 
% ---------------------------------------------------------------------
singleslice         = cfg_menu;
singleslice.tag     = 'singleslice';
singleslice.name    = 'Single Slice Save';
singleslice.val     = {0};
singleslice.help    = {'Choose if you want to save each single slice as matlab structure (you need space on disk for that!)'
    'default: No'};
singleslice.labels = {'No - don''t save each single slice as matlab structure'
               'Yes -  save each single slice as matlab structure'}';
singleslice.values = {0 1};


% ---------------------------------------------------------------------
% raw HARDI Save as matlab structure? 
% ---------------------------------------------------------------------
hardisave         = cfg_menu;
hardisave.tag     = 'hardi';
hardisave.name    = 'HARDI Save';
hardisave.val     = {0};
hardisave.help    = {'Choose if you want to save raw HARDI data as matlab structure (you need space on disk and memory for that!)'
    'default: No'};
hardisave.labels = {'No - don''t save HARDI as matlab structure'
               'Yes -  save HARDI as matlab structure'}';
hardisave.values = {0 1};


% ---------------------------------------------------------------------
% tensor 
% ---------------------------------------------------------------------
tensor         = cfg_exbranch;
tensor.tag     = 'tensor';
tensor.name    = 'Tensor Calculation';
tensor.val     = {filename singleslice hardisave threshold};
tensor.help    = {'Calculates the DTI tensors.'};
tensor.prog    = @dti_tensor_ui;
tensor.vout    = @vout;
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
function dep = vout(job)
dep            = cfg_dep;
dep.sname      = 'DTD data file';
dep.src_output = substruct('.','files');
dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});