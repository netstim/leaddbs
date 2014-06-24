clear
clc
for s=1:7

matlabbatch{1}.spm.util.cat.vols = {
                                    ['/PA/Neuro/_projects/lead/lead/templates/dartel/s',num2str(s+1),'c1dartelmni_6.nii,1']
                                    ['/PA/Neuro/_projects/lead/lead/templates/dartel/s',num2str(s+1),'c2dartelmni_6.nii,1']
                                    ['/PA/Neuro/_projects/lead/lead/templates/dartel/s',num2str(s+1),'c3dartelmni_6.nii,1']
                                    };
matlabbatch{1}.spm.util.cat.name = ['dartelmni_',num2str(6-s),'.nii'];
matlabbatch{1}.spm.util.cat.dtype = 0;

jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear jobs matlabbatch


end