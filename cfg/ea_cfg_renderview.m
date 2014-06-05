% configure file for calculation of tensors
% DTI Configuration file
% BATCH system (Volkmar Glauche)
%
% File created by Andreas Horn


function renderview = ea_cfg_renderview

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
% showatlases Show atlases. 
% ---------------------------------------------------------------------
showatlases         = cfg_menu;
showatlases.tag     = 'showatlases';
showatlases.name    = 'Show Atlases';
showatlases.help    = {'Specify whether to render atlases.'};
showatlases.labels  = {'Yes ? render altases.'
                'No ? do not render altases.'
              }';
showatlases.values  = {1 0};

% ---------------------------------------------------------------------
% showfibres Show fibertracts. 
% ---------------------------------------------------------------------
showfibers         = cfg_menu;
showfibers.tag     = 'showfibers';
showfibers.name    = 'Show Fibers';
showfibers.help    = {'Specify whether to render fibers.'};
showfibers.labels  = {'Yes ? render fibers.'
                'No ? do not render fibers.'
              }';
showfibers.values  = {1 0};

% ---------------------------------------------------------------------
% showconnectivities Show fiber connectivities. 
% ---------------------------------------------------------------------
showconnectivities         = cfg_menu;
showconnectivities.tag     = 'showconnectivities';
showconnectivities.name    = 'Show Fiber connectivities';
showconnectivities.help    = {'Specify whether to render connectivities from electrodes to Atlas regions.'};
showconnectivities.labels  = {'Yes ? render connectivities.'
                'No ? do not render connectivities.'
              }';
showconnectivities.values  = {1 0};

% ---------------------------------------------------------------------
% prolongelectrode visualize longer electrode
% ---------------------------------------------------------------------
prolongelectrode         = cfg_entry;
prolongelectrode.tag     = 'prolongelectrode';
prolongelectrode.name    = 'Prolong electrode rendering';
prolongelectrode.val     = {2};
prolongelectrode.help    = {'Enter the factor to multiply the electrode length with. (Default=2). This is only a visualization option.'};
prolongelectrode.strtype = 'e';
prolongelectrode.num     = [1 1];

% ---------------------------------------------------------------------
% elspec Which Electrode has been used? 
% ---------------------------------------------------------------------
elspec         = cfg_menu;
elspec.tag     = 'elspec';
elspec.name    = 'Electrode Specification';
elspec.help    = {'Select the electrode to render here.'};
elspec.labels  = {'Medtronic 3389'
                'Medtronic 3387'}';
elspec.values  = {'medtronic_3389','medtronic_3387'};

% ---------------------------------------------------------------------
% renderview 
% ---------------------------------------------------------------------
renderview         = cfg_exbranch;
renderview.tag     = 'renderview';
renderview.name    = 'Render View';
renderview.val     = {foldername showatlases showfibers showconnectivities prolongelectrode elspec};
renderview.help    = {'Displays the Electrode reconstruction in 3D.'};
renderview.prog    = @ea_ui_renderview;
renderview.vout    = @vout;
% ---------------------------------------------------------------------
