% Test with fieldtrip headmodel:

% segment volume:
clear
cfg.output='tpm';
  cfg.template       = '/spm8/templates/T1.nii';
  cfg.tpm            = {'/Users/andreashorn/Documents/MATLAB/spm8/tpm/grey.nii','/Users/andreashorn/Documents/MATLAB/spm8/tpm/white.nii','/Users/andreashorn/Documents/MATLAB/spm8/tpm/csf.nii'}
  %cfg.name           = 'tmp';

  % load mri in Fieldtrip-format.
   hdr = spm_vol_nifti('pre_tra.nii');
   img = spm_read_vols(hdr);
   transform = hdr.mat;
   % set up the axes of the volume in voxel coordinates
nx = size(img,1);
ny = size(img,2);
nz = size(img,3);
mri.dim = [nx ny nz];
% store the anatomical data
mri.anatomy = img;
% store the header with all file format specific details
mri.hdr = hdr;

try
  % store the homogenous transformation matrix if present
  mri.transform = transform;
end

try
  % try to determine the units of the coordinate system
  mri = ft_convert_units(mri);
end


seg=ft_volumesegment(cfg,mri);


cfg.method='simbio';


[vol, cfg] = ft_prepare_headmodel(cfg, seg);