ixiids=ea_getIXI_IDs(564);

fsldir = [ea_getearoot, filesep, 'ext_libs', filesep, 'fsl', filesep];
SLICER = ea_getExec([fsldir, 'slicer'], escapePath = 1);


setenv('FSLOUTPUTTYPE','NIFTI');

for ixi=1:length(ixiids)
    [pth, ixibase] = fileparts(ixiids{ixi});
    cmd = [SLICER, ...
           ' ', ixiids{ixi}, filesep, 'glanat.nii.gz' ...
           ' ', ea_space, 't2.nii' ...
           ' -a', ' mono', ixibase, '.png'];

    ea_runcmd(cmd);
end
