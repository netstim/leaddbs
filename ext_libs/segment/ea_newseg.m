function ea_newseg(directory,file,dartel,options)

% we cannot generate the TPM from the SPM TPM anymore since 
% we use the enhanced TPM by Lorio / Draganski:
% http://unil.ch/lren/home/menuinst/data--utilities.html

% if ~exist([options.earoot,'templates',filesep,'TPM.nii'],'file')
%     ea_generate_tpm(options);
% end

disp('Segmentation...');

load([options.earoot,'ext_libs',filesep,'segment',filesep,'segjob12']);
tpminf=[options.earoot,'templates',filesep,'TPM_Lorio_Draganski.nii'];
job.channel.vols{1}=[directory,file,',1'];
for tpm=1:6
    job.tissue(tpm).tpm=[tpminf,',',num2str(tpm)];
    if dartel && tpm < 4
        job.tissue(tpm).native(2)=1;
    end
end
job.warp.write=[1,1];
spm_preproc_run(job); % run "Segment" in SPM 12 (Old "Segment" is now referred to as "Old Segment").


disp('Done.');