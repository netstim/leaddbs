function ea_generate_tpm

spmdir=[fileparts(which('spm')),filesep];
earoot=[fileparts(which('lead')),filesep];

copyfile([spmdir,'toolbox',filesep,'Seg',filesep,'TPM.nii'],[earoot,'templates',filesep,'TPM.nii']);

% spm_file_split([earoot,'templates',filesep,'TPM.nii']);
% 
% 
% 
% %%
% matlabbatch{1}.spm.util.reorient.srcfiles = {[earoot,'templates',filesep,'TPM_00001.nii']
%     [earoot,'templates',filesep,'TPM_00002.nii']
%     [earoot,'templates',filesep,'TPM_00003.nii']
%     [earoot,'templates',filesep,'TPM_00004.nii']
%     [earoot,'templates',filesep,'TPM_00005.nii']
%     [earoot,'templates',filesep,'TPM_00006.nii']
%                                              };
% %%
% matlabbatch{1}.spm.util.reorient.transform.transM = [1.03799350893443 -0.00551554932118825 0.00211932303236285 -0.923782191169304
%                                                      0.00800793946529943 1.0393803580894 -0.00445193813236372 1.02715996003749
%                                                      -0.00568315298822271 -0.0101344712041121 1.01619757821906 0.111266550723032
%                                                      0 0 0 1];
% matlabbatch{1}.spm.util.reorient.prefix = 'w';
% 
% 
% jobs{1}=matlabbatch;
% cfg_util('run',matlabbatch)
% clear jobs matlabbatch
% 
% 
% matlabbatch{1}.spm.util.cat.vols = {[earoot,'templates',filesep,'wTPM_00001.nii,1']
%     [earoot,'templates',filesep,'wTPM_00002.nii,1']
%     [earoot,'templates',filesep,'wTPM_00003.nii,1']
%     [earoot,'templates',filesep,'wTPM_00004.nii,1']
%     [earoot,'templates',filesep,'wTPM_00005.nii,1']
%     [earoot,'templates',filesep,'wTPM_00006.nii,1']
%                                     };
% matlabbatch{1}.spm.util.cat.name = 'TPM.nii';
% matlabbatch{1}.spm.util.cat.dtype = 0;
% 
% jobs{1}=matlabbatch;
% cfg_util('run',matlabbatch)
% clear jobs matlabbatch


%movefile('TPM.nii',[earoot,'templates',filesep,'TPM.nii']);

reslice_nii([earoot,'templates',filesep,'TPM.nii'],[earoot,'templates',filesep,'TPM.nii'],[0.5 0.5 0.5],0,0,3);

% cleanup

delete([earoot,'templates',filesep,'TPM_00001.nii']);
delete([earoot,'templates',filesep,'TPM_00002.nii']);
delete([earoot,'templates',filesep,'TPM_00003.nii']);
delete([earoot,'templates',filesep,'TPM_00004.nii']);
delete([earoot,'templates',filesep,'TPM_00005.nii']);
delete([earoot,'templates',filesep,'TPM_00006.nii']);
delete([earoot,'templates',filesep,'wTPM_00001.nii']);
delete([earoot,'templates',filesep,'wTPM_00002.nii']);
delete([earoot,'templates',filesep,'wTPM_00003.nii']);
delete([earoot,'templates',filesep,'wTPM_00004.nii']);
delete([earoot,'templates',filesep,'wTPM_00005.nii']);
delete([earoot,'templates',filesep,'wTPM_00006.nii']);