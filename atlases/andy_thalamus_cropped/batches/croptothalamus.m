comps={'occipital','postparietal','prefrontal','premotor','primarymotor','sensory','temporal'};
for comp=1:7



matlabbatch{1}.spm.util.imcalc.input = {
                                        '/PA/Neuro/_projects/lead/lead_dbs/atlases/andy_thalamus_cropped/mixed/toccipital.nii,1'
                                        '/PA/Neuro/_projects/lead/lead_dbs/atlases/andy_thalamus_cropped/mixed/tpostparietal.nii,1'
                                        '/PA/Neuro/_projects/lead/lead_dbs/atlases/andy_thalamus_cropped/mixed/tprefrontal.nii,1'
                                        '/PA/Neuro/_projects/lead/lead_dbs/atlases/andy_thalamus_cropped/mixed/tpremotor.nii,1'
                                        '/PA/Neuro/_projects/lead/lead_dbs/atlases/andy_thalamus_cropped/mixed/tprimarymotor.nii,1'
                                        '/PA/Neuro/_projects/lead/lead_dbs/atlases/andy_thalamus_cropped/mixed/tsensory.nii,1'
                                        '/PA/Neuro/_projects/lead/lead_dbs/atlases/andy_thalamus_cropped/mixed/ttemporal.nii,1'
                                        '/PA/Neuro/_projects/lead/lead_dbs/atlases/andy_thalamus_cropped/mixed/thalamus.nii,1'
                                        };
matlabbatch{1}.spm.util.imcalc.output = ['t',comps{comp},'.nii'];
matlabbatch{1}.spm.util.imcalc.outdir = {'/PA/Neuro/_projects/lead/lead_dbs/atlases/andy_thalamus_cropped/mixed/'};
matlabbatch{1}.spm.util.imcalc.expression = ['i',num2str(comp),'.*i8'];
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;



jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear jobs matlabbatch
end