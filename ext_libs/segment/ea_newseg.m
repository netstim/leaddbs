function ea_newseg(directory,file,dartel,options,del)
% SPM NewSegment
% we cannot generate the TPM from the SPM TPM anymore since
% we use the enhanced TPM by Lorio / Draganski:
% http://unil.ch/lren/home/menuinst/data--utilities.html

if nargin < 5
    del = 1;
end

if ~exist([directory, 'c1', file], 'file')
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
    
    % delete unused files
    try delete([directory, 'c4', file]); end
    try delete([directory, 'c5', file]); end
    try delete([directory, 'c6', file]); end
    [~,fn]=fileparts(file);
    try delete([directory,fn,'_seg8.mat']); end
    
    % del==0: keep the deformation field when for ea_normalize_spmnewseg
    if del
        try delete([directory, 'y_', file]); end
        try delete([directory, 'iy_', file]); end
    end
    
    disp('Done.');
else
    disp('Segmentation already done...');
end
