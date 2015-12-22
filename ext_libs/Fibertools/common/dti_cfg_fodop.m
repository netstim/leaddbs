% configure file for operations on streamline trackings results
%
% for BATCH EDITOR system (Volkmar Glauche)
%
% File created by Marco Reisert
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function fodop = dti_cfg_fodop

% ---------------------------------------------------------------------
% filename1 name of HARDI
% ---------------------------------------------------------------------
filename1         = cfg_files;
filename1.tag     = 'filename1';
filename1.name    = 'Select the HARDI-mrstruct';
filename1.help    = {'Select the mrstruct (usually something like blah_HARDI.mat'};
filename1.filter  = 'mat';
filename1.ufilter = '.*';
filename1.num     = [1 1];

% ---------------------------------------------------------------------
% filename1 name of HARDI2
% ---------------------------------------------------------------------
filenameHARDI1         = cfg_files;
filenameHARDI1.tag     = 'filenameHARDI1';
filenameHARDI1.name    = 'Select direction reference as HARDI-mrstruct';
filenameHARDI1.help    = {'Select the mrstruct (usually something like blah_HARDI.mat'};
filenameHARDI1.filter  = 'mat';
filenameHARDI1.ufilter = '.*';
filenameHARDI1.num     = [1 1];


% ---------------------------------------------------------------------
% filename2 name of FOD
% ---------------------------------------------------------------------
filename2         = cfg_files;
filename2.tag     = 'filename2';
filename2.name    = 'Select the FOD-mrstruct';
filename2.help    = {'Select the FOD-mrstruct (usually something like "blah_FOD.mat)"'};
filename2.filter  = 'mat';
filename2.ufilter = '.*';
filename2.num     = [1 1];


% ---------------------------------------------------------------------
% filename2 name of FOD
% ---------------------------------------------------------------------
filename4         = cfg_files;
filename4.tag     = 'filename4';
filename4.name    = 'Select the mrstruct containing the fiber density';
filename4.help    = {'Select the mrstruct (usually something like "blah.mat)"'};
filename4.filter  = 'mat';
filename4.ufilter = '.*';
filename4.num     = [1 1];


% ---------------------------------------------------------------------
% filename2 name of FOD
% ---------------------------------------------------------------------
filenamePIcon         = cfg_files;
filenamePIcon.tag     = 'filenamePIcon';
filenamePIcon.name    = 'Select the PICON file';
filenamePIcon.help    = {'Select the *_PICON.mat file)"'};
filenamePIcon.filter  = 'mat';
filenamePIcon.ufilter = '.*';
filenamePIcon.num     = [1 1];



% ---------------------------------------------------------------------
% filename3 name of nifti
% ---------------------------------------------------------------------
filename3         = cfg_files;
filename3.tag     = 'filename3';
filename3.name    = 'Select a nifti as reference for reslicing';
filename3.help    = {'Select the nifti (usually something like blah.nii)'};
filename3.filter  = 'image';
filename3.ufilter = '.*';
filename3.num     = [1 1];


% ---------------------------------------------------------------------
% filenameROIs name of nifti
% ---------------------------------------------------------------------
filenameROIs         = cfg_files;
filenameROIs.tag     = 'filenameROIs';
filenameROIs.name    = 'Select a nifti containing labels';
filenameROIs.help    = {'Select the nifti (usually something like blah.nii)'};
filenameROIs.filter  = 'image';
filenameROIs.ufilter = '.*';
filenameROIs.num     = [1 1];


% ---------------------------------------------------------------------
% filedef name of deformation field
% ---------------------------------------------------------------------
filedef         = cfg_files;
filedef.tag     = 'filedef';
filedef.name    = 'Select the deformation field';
filedef.help    = {'Select the deformation field (usually created by "New Segment", something like "yblah.nii"'};
filedef.filter  = '.nii';
filedef.ufilter = '.*';
filedef.num     = [1 1];

% ---------------------------------------------------------------------
% filedefi name of inverse deformation field
% ---------------------------------------------------------------------
filedefi         = cfg_files;
filedefi.tag     = 'filedefi';
filedefi.name    = 'Select the inverse deformation field';
filedefi.help    = {'Select the inverse deformation field (usually created by "New Segment", something like "iyblah.nii"'};
filedefi.filter  = '.nii';
filedefi.ufilter = '.*';
filedefi.num     = [1 1];


% ---------------------------------------------------------------------
% thresh, select the threshhold for creating a mask
% ---------------------------------------------------------------------
thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Threshold for creating a mask';
thresh.val     = {0.5};
thresh.help    = {'Threshold for creating a mask of voxels visited by at least this number of fibretracks. Default = 0'};
thresh.strtype = 'e';
thresh.num     = [1 1];



% ---------------------------------------------------------------------
% roiname
% ---------------------------------------------------------------------
roiname         = cfg_files;
roiname.tag     = 'roiname';
roiname.name    = 'White matter mask';
roiname.help    = {'Select a nifti.'};
roiname.filter  = '.mat|.nii';
roiname.ufilter = '.*';
roiname.num     = [1 1];



% ---------------------------------------------------------------------
% Parameters of PItrack
% ---------------------------------------------------------------------
ParametersPI = {'ROI index', ...
                'GMres tolerance', ...
              'GMres maxit', ...
              'GMres restart', ...
              'Local Maxima Threshold', ...
              'Dangular', ...
              'FOD power', ...
              'FOD threshold', ...
              'Length Exponent', ...
              'Boundary Condition',...
              'Normalization Type',...
              'Symmetry', ...
              'Upsampling', ...
              'Speed',
              };
paramtagsPI = {'roisubset', 'tol','maxitgm','restart','thres_toMax','Dang','sh','myeps','kappa','boundary','normtype','sym','upsamp','speed'};
defaultsPI =  {[], 10^-8, 500, 5,0.25, 0.01, 40, 0.01,  0.0 ,1,0,1,[1 1],0};

for k = 1:length(ParametersPI),
    eval([paramtagsPI{k} ' = cfg_entry;']);
    eval([paramtagsPI{k} '.tag = ''' paramtagsPI{k} '''; ']);
    eval([paramtagsPI{k} '.name = ''' ParametersPI{k} '''; ']);
    eval([paramtagsPI{k} '.help = {'' Choose the parameter: ' ParametersPI{k} '''}; ']);
    eval([paramtagsPI{k} '.val = {' num2str(defaultsPI{k}) '};']);
    eval([paramtagsPI{k} '.strtype = ''e'';']);
    eval([paramtagsPI{k} '.num = [1 1];']);
end;
roisubset.num = [0 inf];
roisubset.val = {[]};
upsamp.num = [0 2];
upsamp.val = {[1 1]};

% ---------------------------------------------------------------------
% Parameters of Probtrax
% ---------------------------------------------------------------------
ParametersPX = {'ROI index', ...
                'Walkers per Seed', ...
                'Stepwidth',...
                'Maximal fiber length',...
                'Sigma',...
                'Maximum angle',...
                'Remove length Bias',...
                'No Revisits'};
paramtagsPX = {'roisubset', 'nwalker','step','maxfiblen','sig','maxang','lenbias','norevisits'};
defaultsPX =  {[], 5000,1,1000,0.2,80,1,0};

for k = 1:length(ParametersPX),
    eval([paramtagsPX{k} ' = cfg_entry;']);
    eval([paramtagsPX{k} '.tag = ''' paramtagsPX{k} '''; ']);
    eval([paramtagsPX{k} '.name = ''' ParametersPX{k} '''; ']);
    eval([paramtagsPX{k} '.help = {'' Choose the parameter: ' ParametersPX{k} '''}; ']);
    eval([paramtagsPX{k} '.val = {' num2str(defaultsPX{k}) '};']);
    eval([paramtagsPX{k} '.strtype = ''e'';']);
    eval([paramtagsPX{k} '.num = [1 1];']);
end;
roisubset.num = [0 inf];




roipairs         = cfg_entry;
roipairs.tag     = 'roipairs';
roipairs.name    = 'Select pairs of ROIs to create path trails (a [2 n] array of integers)';
roipairs.val     = {[]};
roipairs.help    = {'some help'};
roipairs.strtype = 'e';
roipairs.num     = [2 inf];



% ---------------------------------------------------------------------
% Parameters of CSD/TFD
% ---------------------------------------------------------------------
Parameters = {'Oversampling Factor for computing TFD', ...
              'Iterations CSD', ...
              'FRF D-axial', ...
              'FRF D-radial', ...
              'L1-Regularization Strength',...
              'Smoothing power',...
              'Number Directions',...
              'Tensororder TFD',...
              'Global Tikhohnov TFD',...
              'Boundary Variation Penalty'};
paramtags = {'osamp','maxit','Dax','Drad','lambdaL1','powsm','numdir','planorder','lambda1','lambda2' };
defaults =  {2,      200, 1, 0.15,   50,      14,     128,      1,         0, 10^(-4) };

for k = 1:length(Parameters),
    eval([paramtags{k} ' = cfg_entry;']);
    eval([paramtags{k} '.tag = ''' paramtags{k} '''; ']);
    eval([paramtags{k} '.name = ''' Parameters{k} '''; ']);
    eval([paramtags{k} '.help = {'' Choose the parameter: ' Parameters{k} '''}; ']);
    eval([paramtags{k} '.val = {' num2str(defaults{k}) '};']);
    eval([paramtags{k} '.strtype = ''e'';']);
    eval([paramtags{k} '.num = [1 1];']);
end;
 



% ---------------------------------------------------------------------
% Parameters of CSA
% ---------------------------------------------------------------------
ParametersCSA = {'SH cutoff', ...
              'Laplace-Beltrami penalty', ...
              'Number of Directions (output space)',...
              'FRF D-axial'};
paramtagsCSA = {'shcutoff','laplacebeltrami','numdirsout','Dax'};
defaultsCSA =  {6, 0.006,  128, 1 };

for k = 1:length(ParametersCSA),
    eval([paramtagsCSA{k} ' = cfg_entry;']);
    eval([paramtagsCSA{k} '.tag = ''' paramtagsCSA{k} '''; ']);
    eval([paramtagsCSA{k} '.name = ''' ParametersCSA{k} '''; ']);
    eval([paramtagsCSA{k} '.help = {'' Choose the parameter: ' ParametersCSA{k} '''}; ']);
    eval([paramtagsCSA{k} '.val = {' num2str(defaultsCSA{k}) '};']);
    eval([paramtagsCSA{k} '.strtype = ''e'';']);
    eval([paramtagsCSA{k} '.num = [1 1];']);
end;
 
ParametersCSDFC = {'Number Iterations', ...
              'Number CG rounds', ...
              'Number of Directions (output space)',...
              'FRF D-axial',...
              'lambda AFC',...
              'alpha curvature',...
              'lambda FC (symmetric)',...
              'lambda FC (asymmetric)',...
              'lambda isotropic (symmetric)',...
              'lambda isotropic (asymmetric)',...
              'lambda Laplace Beltrami',...
              'lambda Tikhonov (symmetric)',...
              'lambda Tikhonov (asymmetric)',...
              'symmetric output',...
              };
paramtagsCSDFC = {'numit','numround','numoutdir','Daxial','lambda','alpha','purefc_sym','purefc_asym',...
                 'iso_sym','iso_asym','lambda_lb','gamma_sym','gamma_asym','symout'};
defaultsCSDFC =  {50,4,128,1,0,1,0.005,0,0,0,0.0001,0,0,1};
defaultsCSDAFC =  {50,4,128,1,0.005,1,0,0,0,0,0.0001,0,0,1};

for k = 1:length(ParametersCSDFC),
    eval([paramtagsCSDFC{k} ' = cfg_entry;']);
    eval([paramtagsCSDFC{k} '.tag = ''' paramtagsCSDFC{k} '''; ']);
    eval([paramtagsCSDFC{k} '.name = ''' ParametersCSDFC{k} '''; ']);
    eval([paramtagsCSDFC{k} '.help = {'' Choose the parameter: ' ParametersCSDFC{k} '''}; ']);
    eval([paramtagsCSDFC{k} '.val = {' num2str(defaultsCSDFC{k}) '};']);
    eval([paramtagsCSDFC{k} '.strtype = ''e'';']);
    eval([paramtagsCSDFC{k} '.num = [1 1];']);
end;






% ---------------------------------------------------------------------
% roidef 
% ---------------------------------------------------------------------
roidef         = cfg_branch;
roidef.tag     = 'roidef';
roidef.name    = 'Definition of area of Reconstruction';
roidef.help    = {'Choose a coregistered Nifti containing the probabilistiv white matter segmentation.'};
roidef.val     = {roiname thresh};

% ---------------------------------------------------------------------
% newfilename name of resulting file
% ---------------------------------------------------------------------
newfilename         = cfg_entry;
newfilename.tag     = 'newfilename';
newfilename.name    = 'New File Name (mrStruct)';
newfilename.val     = {'.mat'};
newfilename.help    = {'Type in the name of the new file, if you leave this empty a default file name is generated'};
newfilename.strtype = 's';
newfilename.num     = [1 Inf];

newfilename2         = cfg_entry;
newfilename2.tag     = 'newfilename';
newfilename2.name    = 'New File Name (Nifti)';
newfilename2.val     = {'.nii'};
newfilename2.help    = {'Type in the name of the new file, if you leave this empty a default file name is generated'};
newfilename2.strtype = 's';
newfilename2.num     = [1 Inf];

newfilenamemat         = cfg_entry;
newfilenamemat.tag     = 'newfilenamemat';
newfilenamemat.name    = 'New File Name (.mat)';
newfilenamemat.val     = {'.mat'};
newfilenamemat.help    = {'Type in the name of the file contaning connectivity information. If you leave this empty a default file name is generated'};
newfilenamemat.strtype = 's';
newfilenamemat.num     = [1 Inf];


interp         = cfg_menu;
interp.tag     = 'interp';
interp.name    = 'Interpolation';
interp.help    = {'The method by which the images are sampled when being written in a different space. Nearest Neighbour is fastest, but not normally recommended. It can be useful for re-orienting images while preserving the original intensities (e.g. an image consisting of labels). Bilinear Interpolation is OK for PET, or realigned and re-sliced fMRI. If subject movement (from an fMRI time series) is included in the transformations then it may be better to use a higher degree approach. Note that higher degree B-spline interpolation/* \cite{thevenaz00a,unser93a,unser93b}*/ is slower because it uses more neighbours.'};
interp.labels = {
                 'Nearest neighbour'
                 'Trilinear'
                 '2nd Degree B-Spline'
                 '3rd Degree B-Spline'
                 '4th Degree B-Spline'
                 '5th Degree B-Spline'
                 '6th Degree B-Spline'
                 '7th Degree B-Spline'
}';
interp.values = {0 1 2 3 4 5 6 7};
interp.def = @(x) 1;





sitype         = cfg_menu;
sitype.tag     = 'sitype';
sitype.name    = 'Type of scalar map';
sitype.help    = {'Type of scalar maps which is computed from the FOD.'};
sitype.labels = {
                 'Generalized Fractional Anisotropy'
                 'mean FOD'                 
                 'mean FOD (b0 normalized)'
                 'maximum FOD'
                 'maximum FOD (b0 normalized)'                 
}';
sitype.values = {0 1 2 3 4};
sitype.def = @(x) 0;




modifyFOD         = cfg_menu;
modifyFOD.tag     = 'modifyFOD';
modifyFOD.name    = 'Scale FOD by TFD values';
modifyFOD.help    = {'The original FOD scaling is overwritten by TFD.'};
modifyFOD.labels = {
                 'Yes'
                 'No'
}';
modifyFOD.values = {1 0};
modifyFOD.def = @(x) 0;



savePImaps         = cfg_menu;
savePImaps.tag     = 'savePImaps';
savePImaps.name    = 'Save PImaps to generate Path trails.';
savePImaps.help    = {'Information about path trails is saves'};
savePImaps.labels = {
                 'Yes'
                 'No'
}';
savePImaps.values = {1 0};
savePImaps.def = @(x) 0;










noutdir         = cfg_menu;
noutdir.tag     = 'noutdir';
noutdir.name    = 'Number of Output Directions';
noutdir.help    = {'Number of directions in output space'};
noutdir.labels = {
                 '32'
                 '64'                 
                 '128'
}';
noutdir.values = {32 64 128};
noutdir.def = @(x) 64;



modulate       = cfg_menu;
modulate.tag     = 'modulate';
modulate.name    = 'Modulation';
modulate.help    = {'The method with which the values are scaled. To preserve surface integrals a FOD is scaled by |Jn|^4 /det(J)^2, where J is the Jacobian of the deformation field. To preserve fiber density hte FOD is scaled by |Jn|^3 /det(J).'};
modulate.labels = {
                 'No modulation'
                 'Surface integral preserving modulation (Raeffel 2012)'
                 'Fiber Density Preserving modulation (Hong 2009)'
}';
modulate.values = {0 1 2};
modulate.def = @(x) 1;


axialsigma         = cfg_entry;
axialsigma.tag     = 'axialsigma';
axialsigma.name    = 'Gaussian width in axial direction (in voxel units)';
axialsigma.val     = {3};
axialsigma.help    = {'Smoothing width in axial direction'};
axialsigma.strtype = 'e';
axialsigma.num     = [1 1];


radialsigma         = cfg_entry;
radialsigma.tag     = 'radialsigma';
radialsigma.name    = 'Gaussian width in radial direction (in voxel units)';
radialsigma.val     = {1};
radialsigma.help    = {'Smoothing width in radial direction'};
radialsigma.strtype = 'e';
radialsigma.num     = [1 1];

noiselevel         = cfg_entry;
noiselevel.tag     = 'noiselevel';
noiselevel.name    = 'Noise level (in signal units)';
noiselevel.val     = {20};
noiselevel.help    = {'Noise level used for MLE estimate (0 means level is estimated)'};
noiselevel.strtype = 'e';
noiselevel.num     = [1 1];

ksz         = cfg_entry;
ksz.tag     = 'ksz';
ksz.name    = 'Radius of kernel in voxel units';
ksz.val     = {2};
ksz.help    = {'size of filter cube'};
ksz.strtype = 'e';
ksz.num     = [1 1];



% ---------------------------------------------------------------------
% fiout (nifti or mrstruct) with reslicing option
% ---------------------------------------------------------------------


nifti         = cfg_exbranch;
nifti.tag     = 'newfinamenifti';
nifti.name    = 'New File Name (Nifti)';
nifti.help    = {'You can choose a reference image for reslicing'};
nifti.val     = {newfilename2};


niftireslice         = cfg_exbranch;
niftireslice.tag     = 'newfinameniftireslice';
niftireslice.name    = 'New File Name (Nifti with Reslicing)';
niftireslice.help    = {'Choose a filename'};
niftireslice.val     = {filename3 newfilename2};



fiout         = cfg_choice;
fiout.tag     = 'fiout';
fiout.name    = 'Output File';
fiout.help    = {'Select mrStruct or Nifti as output type'};
fiout.values  = {newfilename nifti niftireslice};



% ---------------------------------------------------------------------
% createFOD by CSD-L1
% ---------------------------------------------------------------------

newfilenameFOD = newfilename; newfilenameFOD.name  = 'New File Name FOD (mrStruct)';
createbyCSD         = cfg_exbranch;
createbyCSD.tag     = 'createbyCSD';
createbyCSD.name    = 'Create FOD by L1-based CSD';
createbyCSD.help    = {'Create Fiber Orientation Distribution with CSD.'};
pp  = (cellfun(@(x) [x ' '],paramtags,'uniformoutput',false)); pp = cat(2,pp{2:7});
createbyCSD.val     = eval(['{filename1 roidef newfilenameFOD ' pp ' }']);
createbyCSD.prog    = @(job)dti_fodop_ui('doL1CSD',job);
createbyCSD.vout    = @vout;



% ---------------------------------------------------------------------
% createFOD by CSA (Aganj2009)
% ---------------------------------------------------------------------

newfilenameFOD = newfilename; newfilenameFOD.name  = 'New File Name FOD (mrStruct)';
createbyCSA         = cfg_exbranch;
createbyCSA.tag     = 'createbyCSA';
createbyCSA.name    = 'Create FOD by Constant Solid Angle (CSA) approach';
createbyCSA.help    = {'Create Fiber Orientation Distribution with CSA approach (Aganj2009).'};
pp  = (cellfun(@(x) [x ' '],paramtagsCSA,'uniformoutput',false)); pp = cat(2,pp{1:3});
createbyCSA.val     = eval(['{filename1 newfilenameFOD ' pp ' }']);
createbyCSA.prog    = @(job)dti_fodop_ui('doCSA',job);
createbyCSA.vout    = @vout;


% ---------------------------------------------------------------------
% createFOD by CSD (Tournier)
% ---------------------------------------------------------------------

newfilenameFOD = newfilename; newfilenameFOD.name  = 'New File Name FOD (mrStruct)';
createbyCSDt         = cfg_exbranch;
createbyCSDt.tag     = 'createbyCSDTournier';
createbyCSDt.name    = 'Create FOD by CSD (Tournier)';
createbyCSDt.help    = {'Create Fiber Orientation Distribution with CSD approach (Tournier).'};
pp  = (cellfun(@(x) [x ' '],paramtagsCSA,'uniformoutput',false)); pp = cat(2,pp{[1 3 4]});
createbyCSDt.val     = eval(['{filename1 newfilenameFOD roidef ' pp ' }']);
createbyCSDt.prog    = @(job)dti_fodop_ui('doCSDTournier',job);
createbyCSDt.vout    = @vout;



% ---------------------------------------------------------------------
% createFOD by CSD with fiber continuity
% ---------------------------------------------------------------------


createbyCSDFC         = cfg_exbranch;
createbyCSDFC.tag     = 'createbyCSFC';
createbyCSDFC.name    = 'Create FOD by CSD-Fiber Continuity';
createbyCSDFC.help    = {'Create Fiber Orientation Distribution by CSD with Fiber Continuity regularizer (Reisert).'};
pp  = (cellfun(@(x) [x ' '],paramtagsCSDFC,'uniformoutput',false)); pp = cat(2,pp{:});
createbyCSDFC.val     = eval(['{filename1 newfilenameFOD roidef ' pp ' }']);
createbyCSDFC.prog    = @(job)dti_fodop_ui('doCSDFC',job);
createbyCSDFC.vout    = @vout;





% ---------------------------------------------------------------------
% createFODTFD
% ---------------------------------------------------------------------
fioutTFD = fiout; fioutTFD.name = 'New File Name TFD (mrStruct or Nifti)';
createbyCSDTFD         = cfg_exbranch;
createbyCSDTFD.tag     = 'createbyCSDTFD';
createbyCSDTFD.name    = 'Create FOD by L1-CSD with Tensor Divergence';
createbyCSDTFD.help    = {'Create Fiber Orientation Distribution with CSD, where mean density per voxel is determined by  TFD.'};
pp  = (cellfun(@(x) [x ' '],paramtags,'uniformoutput',false)); pp = cat(2,pp{:});
createbyCSDTFD.val     = eval(['{filename1 roidef newfilenameFOD fioutTFD ' pp ' }']);
createbyCSDTFD.prog    = @(job)dti_fodop_ui('doL1CSD_TFD',job);
createbyCSDTFD.vout    = @voutFODTFD;


% ---------------------------------------------------------------------
% createTFD
% ---------------------------------------------------------------------
filename1FOD = filename1;
filename1FOD.name    = 'Select the FOD-mrstruct';
filename1FOD.help    = {'Select the mrstruct (usually something like blah_FOD.mat'};
createbyTFD         = cfg_exbranch;
createbyTFD.tag     = 'createbyCSDTFD';
createbyTFD.name    = 'Create TFD from FOD';
createbyTFD.help    = {'Create tensor fiber density from given FOD.'};
pp  = (cellfun(@(x) [x ' '],paramtags,'uniformoutput',false)); pp = cat(2,pp{[1 8:10]});
createbyTFD.val     = eval(['{filename1FOD roidef  fioutTFD modifyFOD ' pp ' }']);
createbyTFD.prog    = @(job)dti_fodop_ui('doTFD',job);
createbyTFD.vout    = @voutTFD;


% ---------------------------------------------------------------------
% deformFOD
% ---------------------------------------------------------------------
deformFOD         = cfg_exbranch;
deformFOD.tag     = 'deformFOD';
deformFOD.name    = 'Deformation of FOD';
deformFOD.help    = {'Non Rigid deformation of FODs by deformation fields created with NewSegment.'};
deformFOD.val     = {filename2 filedef filedefi interp modulate noutdir newfilename};
deformFOD.prog    = @(job)dti_fodop_ui('deformFOD',job);
deformFOD.vout    = @vout;



% ---------------------------------------------------------------------
% deformFOD
% ---------------------------------------------------------------------
deformFDalong         = cfg_exbranch;
deformFDalong.tag     = 'deformFDalong';
deformFDalong.name    = 'Deformation of FD along FOD';
deformFDalong.help    = {'Non Rigid deformation of FD intepreted as mean of FOD.'};
deformFDalong.val     = {filename2 filename4 filedef filedefi interp modulate newfilename};
deformFDalong.prog    = @(job)dti_fodop_ui('deformFDalong',job);
deformFDalong.vout    = @vout;


% ---------------------------------------------------------------------
% smoothFOD
% ---------------------------------------------------------------------
smoothFOD         = cfg_exbranch;
smoothFOD.tag     = 'smoothFOD';
smoothFOD.name    = 'Anisotropic Smooth of FOD/HARDI';
smoothFOD.help    = {'Smooth along the current FOD direction.'};
smoothFOD.val     = {filename2  axialsigma radialsigma  newfilename};
smoothFOD.prog    = @(job)dti_fodop_ui('smoothFOD',job);
smoothFOD.vout    = @vout;


% ---------------------------------------------------------------------
% smoothFOD
% ---------------------------------------------------------------------
smoothFODRic         = cfg_exbranch;
smoothFODRic.tag     = 'smoothFODRic';
smoothFODRic.name    = 'Anisotropic Smooth of FOD/HARDI (Rician)';
smoothFODRic.help    = {'Smooth along the current FOD direction.'};
smoothFODRic.val     = {filename2  axialsigma radialsigma ksz noiselevel newfilename};
smoothFODRic.prog    = @(job)dti_fodop_ui('smoothFODRic',job);
smoothFODRic.vout    = @vout;



% ---------------------------------------------------------------------
% get mFOD
% ---------------------------------------------------------------------
getsiFOD         = cfg_exbranch;
getsiFOD.tag     = 'getsiFOD';
getsiFOD.name    = 'Compute scalar map from FOD';
getsiFOD.help    = {'Computes several types of scalar indices from FOD.'};
getsiFOD.val     = {filename2 sitype fiout};
getsiFOD.prog    = @(job)dti_fodop_ui('getsiFOD',job);
getsiFOD.vout    = @voutSI;



% ---------------------------------------------------------------------
% resample HARDI
% ---------------------------------------------------------------------

resampleHARDI         = cfg_exbranch;
resampleHARDI.tag     = 'resampleHARDI';
resampleHARDI.name    = 'Resample Diffusion directions of HARDI data';
resampleHARDI.help    = {'Resamples Orientation Information.'};
resampleHARDI.val     = {filename1 filenameHARDI1};
resampleHARDI.prog    = @(job)dti_fodop_ui('resampleHARDI',job);
resampleHARDI.vout    = @voutHARDI;


% ---------------------------------------------------------------------
% pathIntegral tracking
% ---------------------------------------------------------------------






outftr         = cfg_branch;
outftr.tag     = 'outftr';
outftr.name    = 'Specify name of FTR';
outftr.val     = {newfilenamemat};
outftr.help    = {'Specify a output filename.'};

% ---------------------------------------------------------------------
% autoimg Automatically generate Output Filename
% ---------------------------------------------------------------------
noftr        = cfg_const;
noftr.tag    = 'noftr';
noftr.name   = 'Produce no fiber tracks';
noftr.val    = {true};
noftr.help   = {'Produce no fiber tracks by streamlining.'};

% ---------------------------------------------------------------------
% outname Output Naming Scheme
% ---------------------------------------------------------------------
streamftr         = cfg_choice;
streamftr.tag     = 'streamftr';
streamftr.name    = 'Produce Streamlines (single ROI)';
streamftr.values  = {outftr noftr};
streamftr.help    = {'Produce Streamlines'};

streamftrpaired         = cfg_choice;
streamftrpaired.tag     = 'streamftrpaired';
streamftrpaired.name    = 'Produce Streamlines (ROI pairs)';
streamftrpaired.values  = {outftr noftr};
streamftrpaired.help    = {'Produce Streamlines.'};




PItrack         = cfg_exbranch;
PItrack.tag     = 'PItrack';
PItrack.name    = 'Path Integral Connectivity';
PItrack.help    = {'please help me out.'};
pp  = (cellfun(@(x) [x ' '],paramtagsPI,'uniformoutput',false)); pp = cat(2,pp{:});
PItrack.val     = eval(['{filename2 roiname thresh filenameROIs newfilenamemat savePImaps streamftr streamftrpaired ' pp ' }']);
PItrack.prog    = @(job)dti_fodop_ui('PItrack',job);
PItrack.vout    = @voutPITRACK;


PItrails         = cfg_exbranch;
PItrails.tag     = 'PItrails';
PItrails.name    = 'Create PI trails';
PItrails.val     = {filenamePIcon roipairs};
PItrails.prog    = @(job)dti_fodop_ui('PItrails',job);
PItrails.vout    = @voutPITRAILS;



% ---------------------------------------------------------------------
% probtrax
% ---------------------------------------------------------------------

PXtrack         = cfg_exbranch;
PXtrack.tag     = 'PXtrack';
PXtrack.name    = 'Probtrackx';
PXtrack.help    = {'please help me out.'};
pp  = (cellfun(@(x) [x ' '],paramtagsPX,'uniformoutput',false)); pp = cat(2,pp{:});
PXtrack.val     = eval(['{filename2 roiname thresh filenameROIs newfilenamemat ' pp ' }']);
PXtrack.prog    = @(job)dti_fodop_ui('PXtrack',job);
PXtrack.vout    = @voutPXTRACK;



% ---------------------------------------------------------------------
% mapop 
% ---------------------------------------------------------------------
fodop         = cfg_choice;
fodop.tag     = 'fodop';
fodop.name    = 'Operations with FODs';
fodop.help    = {'Operations with fiber orientation distributions (FOD).'};
fodop.values  = {createbyCSD createbyCSDt createbyCSDFC createbyCSA createbyTFD deformFOD smoothFOD smoothFODRic getsiFOD resampleHARDI PItrack PItrails PXtrack};

% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
function dep = vout(job)
dep            = cfg_dep;
dep.sname      = 'FOD mrstruct';
dep.src_output = substruct('.','FOD');
dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});

function dep = voutPITRACK(job)
dep            = cfg_dep;
dep.sname      = 'PI Connectivity Information and (optionally) Maps';
dep.src_output = substruct('.','PIinfo');
dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});

function dep = voutPXTRACK(job)
dep            = cfg_dep;
dep.sname      = 'ProbtrackX Connectivity Information';
dep.src_output = substruct('.','PXinfo');
dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});


function dep = voutPITRAILS(job)
dep            = cfg_dep;
dep.sname      = 'Nifti with PI-trail maps';
dep.src_output = substruct('.','mapnames');
dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});


function dep = voutHARDI(job)
dep            = cfg_dep;
dep.sname      = 'HARDI mrstruct';
dep.src_output = substruct('.','HARDI');
dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});


function dep = voutSI(job)
dep            = cfg_dep;
dep.sname      = 'Scalar Index (nifti/mrstruct)';
dep.src_output = substruct('.','SIndex');
dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});


function dep = voutFODTFD(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'FOD mrstruct';
dep(1).src_output = substruct('.','FOD');
dep(1).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'TFD (nifti/mrstruct)';
dep(2).src_output = substruct('.','TFD');
dep(2).tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});


function dep = voutTFD(job)
dep            = cfg_dep;
dep.sname      = 'TFD mrstruct';
dep.src_output = substruct('.','TFD');
dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});
