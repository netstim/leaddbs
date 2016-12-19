% configure file for operations on streamline trackings results
%
% for BATCH EDITOR system (Volkmar Glauche)
%
% File created by Susanne Schnell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function roiop = dti_cfg_roiop

% ---------------------------------------------------------------------
% newfilename name of resulting mask after operation
% ---------------------------------------------------------------------
newfilename         = cfg_entry;
newfilename.tag     = 'newfilename';
newfilename.name    = 'New File Name';
newfilename.val     = {'.mat'};
newfilename.help    = {'Type in the name of the new mask file, if you leave this empty a default file name is generated.'};
newfilename.strtype = 's';
newfilename.num     = [1 Inf];

% ---------------------------------------------------------------------
% coordinated of the starting point for spherical growing
% ---------------------------------------------------------------------
coords         = cfg_entry;
coords.tag     = 'coords';
coords.name    = 'coords';
coords.val     = {[64 64 10]};
coords.help    = {'Enter the coordinates for the midpoint of the sphere (x y z).'};
coords.strtype = 'e';
coords.num     = [1 3];

% ---------------------------------------------------------------------
% number of voxels
% ---------------------------------------------------------------------
sizeSphere         = cfg_entry;
sizeSphere.tag     = 'sizeSphere';
sizeSphere.name    = 'sizeSphere';
sizeSphere.val     = {1};
sizeSphere.help    = {'Enter the wanted size of the sphere to grow.'};
sizeSphere.strtype = 'e';
sizeSphere.num     = [1 1];

% ---------------------------------------------------------------------
% number of voxels
% ---------------------------------------------------------------------
voxels         = cfg_entry;
voxels.tag     = 'voxels';
voxels.name    = 'voxels';
voxels.val     = {1};
voxels.help    = {'Enter the number of voxels for Erosion or Dilation operation.'};
voxels.strtype = 'e';
voxels.num     = [1 1];

% ---------------------------------------------------------------------
% mask, select the correct mask by a number
% ---------------------------------------------------------------------
mask2         = cfg_entry;
mask2.tag     = 'mask2';
mask2.name    = 'second mask number';
mask2.help    = {'Select the mask by a number. You need to remember the mask position by yourself.'};
mask2.strtype = 'e';
mask2.num     = [1 1];

% ---------------------------------------------------------------------
% mask, select the correct mask by a number
% ---------------------------------------------------------------------
mask1         = cfg_entry;
mask1.tag     = 'mask1';
mask1.name    = 'mask number';
mask1.help    = {'Select the mask by a number. You need to remember the mask position by yourself.'};
mask1.strtype = 'e';
mask1.num     = [1 1];

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
% maskname name of resulting operation
% ---------------------------------------------------------------------
maskname         = cfg_entry;
maskname.tag     = 'maskname';
maskname.name    = 'New Mask Name';
maskname.val     = {'new_ROI'};
maskname.help    = {'Type in the name of the new ROI, if you leave this empty a default file name is generated.'};
maskname.strtype = 's';
maskname.num     = [1 Inf];

% ---------------------------------------------------------------------
% GrowSphere in ROI from point, input of maskstruct, masknumber, size in mm, the new maskname and the new filename
% ---------------------------------------------------------------------
GrowSphere         = cfg_exbranch;
GrowSphere.tag     = 'GrowSphere';
GrowSphere.name    = 'GrowSphere in ROI';
GrowSphere.help    = {'Let grow a sphere from a defined point, give the size in mm.'};
GrowSphere.val     = {roiname maskname sizeSphere coords newfilename};
GrowSphere.prog    = @(job)dti_roiop_ui('GrowSphere',job);
GrowSphere.vout    = @vout;

% ---------------------------------------------------------------------
% Dilation of ROI, input of maskstruct, masknumber, number of voxels, the new maskname and the new filename
% ---------------------------------------------------------------------
Dilation         = cfg_exbranch;
Dilation.tag     = 'Dilation';
Dilation.name    = 'Dilation of ROI';
Dilation.help    = {'Make a dilation operation on the selected mask, give the number of voxels to dilate. This is a standard morphological operation.'};
Dilation.val     = {roiname mask1 maskname voxels newfilename};
Dilation.prog    = @(job)dti_roiop_ui('Dilation',job);
Dilation.vout    = @vout;

% ---------------------------------------------------------------------
% Erosion of ROI, input of maskstruct, masknumber, number of voxels, the new maskname and the new filename
% ---------------------------------------------------------------------
Erosion         = cfg_exbranch;
Erosion.tag     = 'Erosion';
Erosion.name    = 'Erosion of ROI';
Erosion.help    = {'Make an Erosion operation on the selected mask, give the number of voxels to erose. This is a standard morphological operation.'};
Erosion.val     = {roiname mask1 maskname voxels newfilename};
Erosion.prog    = @(job)dti_roiop_ui('Erosion',job);
Erosion.vout    = @vout;

% ---------------------------------------------------------------------
% XORop input of two ROIs, masknumbers, the new maskname and the new filename
% ---------------------------------------------------------------------
XORop         = cfg_exbranch;
XORop.tag     = 'XORop';
XORop.name    = 'XOR operation on two masks';
XORop.help    = {'Select two mask which you want to combine with logical XOR operation.'};
XORop.val     = {roiname mask1 mask2 maskname newfilename};
XORop.prog    = @(job)dti_roiop_ui('XORop',job);
XORop.vout    = @vout;

% ---------------------------------------------------------------------
% ORop input of two ROIs, masknumbers, the new maskname and the new filename
% ---------------------------------------------------------------------
ORop         = cfg_exbranch;
ORop.tag     = 'ORop';
ORop.name    = 'OR operation on two masks';
ORop.help    = {'Select two mask which you want to combine with logical OR operation.'};
ORop.val     = {roiname mask1 mask2 maskname newfilename};
ORop.prog    = @(job)dti_roiop_ui('ORop',job);
ORop.vout    = @vout;

% ---------------------------------------------------------------------
% ANDop input of two ROIs, masknumbers, the new maskname and the new filename
% ---------------------------------------------------------------------
ANDop         = cfg_exbranch;
ANDop.tag     = 'ANDop';
ANDop.name    = 'AND operation on two masks';
ANDop.help    = {'Select two mask which you want to combine with logical AND operation.'};
ANDop.val     = {roiname mask1 mask2 maskname newfilename};
ANDop.prog    = @(job)dti_roiop_ui('ANDop',job);
ANDop.vout    = @vout;

% ---------------------------------------------------------------------
% Invert ROI, input of maskstruct, masknumber, new maskname and the new filename
% ---------------------------------------------------------------------
Invert        = cfg_exbranch;
Invert.tag     = 'Invert';
Invert.name    = 'Invert ROI';
Invert.help    = {'Invert the mask (1 becomes 0 and 0 becomes 1).'};
Invert.val     = {roiname mask1 maskname newfilename};
Invert.prog    = @(job)dti_roiop_ui('Invert',job);
Invert.vout    = @vout;

% ---------------------------------------------------------------------
% roiop 
% ---------------------------------------------------------------------
roiop         = cfg_choice;
roiop.tag     = 'roiop';
roiop.name    = 'Operations with ROIs (maskstructs)';
roiop.help    = {'Operations with ROIs (maskstructs).'};
roiop.values  = {Invert ANDop ORop XORop Erosion Dilation GrowSphere};

% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
function dep = vout(job)
dep            = cfg_dep;
dep.sname      = 'MASKstruct After Operation';
dep.src_output = substruct('.','files');
dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});