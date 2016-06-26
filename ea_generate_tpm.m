function ea_generate_tpm(varargin)

spmdir=[fileparts(which('spm')),filesep];
earoot=[ea_getearoot];
if nargin
    earoot=varargin{1}.earoot;
end

switch spm('ver')
    case 'SPM8'
        copyfile([spmdir,'toolbox',filesep,'Seg',filesep,'TPM.nii'],[earoot,'templates',filesep,'TPM.nii']);
    case 'SPM12'
        copyfile([spmdir,'tpm',filesep,'TPM.nii'],[earoot,'templates',filesep,'TPM.nii']);
end


% for tp=1:6
%     matlabbatch{1}.spm.util.imcalc.input = {[earoot,'templates',filesep,'mni_hires.nii']
%         [earoot,'templates',filesep,'TPM.nii,1']
%         [earoot,'templates',filesep,'TPM.nii,2']
%         [earoot,'templates',filesep,'TPM.nii,3']
%         [earoot,'templates',filesep,'TPM.nii,4']
%         [earoot,'templates',filesep,'TPM.nii,5']
%         [earoot,'templates',filesep,'TPM.nii,6']};
%     matlabbatch{1}.spm.util.imcalc.output = [earoot,'templates',filesep,'TPM_',num2str(tp),'.nii'];
%     matlabbatch{1}.spm.util.imcalc.outdir = {''};
%     matlabbatch{1}.spm.util.imcalc.expression = ['i',num2str(tp+1)];
%     matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
%     matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
%     matlabbatch{1}.spm.util.imcalc.options.mask = 0;
%     matlabbatch{1}.spm.util.imcalc.options.interp = 1;
%     matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
%     cfg_util('run',{matlabbatch});
%     fis{tp}=[earoot,'templates',filesep,'TPM_',num2str(tp),'.nii'];
% end
% clear matlabbatch
% matlabbatch{1}.spm.util.cat.vols = fis;
% matlabbatch{1}.spm.util.cat.name = [earoot,'templates',filesep,'TPM.nii'];
% matlabbatch{1}.spm.util.cat.dtype = 4;
% cfg_util('run',{matlabbatch});
%
% for fi=1:length(fis)
%     delete(fis{fi});
% end
% delete([earoot,'templates',filesep,'TPM.mat']);
