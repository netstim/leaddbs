ixiids=ea_getIXI_IDs(564);

fsldir = [ea_getearoot, filesep, 'ext_libs', filesep, 'fsl', filesep];
if ispc
    SLICER = [fsldir, 'slicer.exe'];
else 
    SLICER = [fsldir, 'slicer.', computer('arch')];    
end

setenv('FSLOUTPUTTYPE','NIFTI');

for ixi=1:length(ixiids)
    [pth, ixibase] = fileparts(ixiids{ixi});
    cmd = [SLICER, ...
           ' ', ixiids{ixi}, filesep, 'glanat.nii.gz' ...
           ' ', ea_space, 't2.nii' ...
           ' -a', ' mono', ixibase, '.png'];
    
    if ~ispc
        system(['bash -c "', cmd, '"']);
    else
        system(cmd);    
    end
end
