ixiids=ea_getIXI_IDs(564);
setenv('FSLDIR','/usr/local/fsl')
!/usr/local/fsl/etc/fslconf/fsl.sh


    % set FSL environment
setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be

for ixi=1:length(ixiids)
[pth,ixibase]=fileparts(ixiids{ixi});
system(['/usr/local/fsl/bin/slicer ',ixiids{ixi},filesep,'glanat.nii.gz /PA/Neuro/_projects/lead/lead_dbs/templates/mni_hires.nii -a mono',ixibase,'.png']);
end