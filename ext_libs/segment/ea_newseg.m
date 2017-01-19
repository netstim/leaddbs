function ea_newseg(directory,files,dartel,options,del,force)
% SPM NewSegment
% we cannot generate the TPM from the SPM TPM anymore since
% we use the enhanced TPM by Lorio / Draganski:
% http://unil.ch/lren/home/menuinst/data--utilities.html

% keep the 'y_*' and 'iy_*' file by default
if ~exist('del','var')
    del = 0;
end
if ~exist('force','var')
    force=0;
end
if ischar(files)
    files={files};
end

if ea_checktpmresolution
    ea_create_tpm_darteltemplate;
end

if (~dartel && exist([directory, 'c1', files{1}], 'file') || dartel && exist([directory, 'rc1', files{1}], 'file')) && (~force)
    disp('Segmentation already done!');
else
    disp('Segmentation...');

    load([options.earoot,'ext_libs',filesep,'segment',filesep,'segjob12']);
    tpminf=[ea_space,'TPM.nii'];
    for fi=1:length(files) % add channels for multispectral segmentation
        job.channel(fi).vols = {[directory,files{fi},',1']};
        job.channel(fi).biasreg = 0.001;
        job.channel(fi).biasfwhm = 60;
        job.channel(fi).write = [0 0];
    end
    tn=ea_open_vol(tpminf);
    for tpm=1:length(tn)
        job.tissue(tpm).tpm=[tpminf,',',num2str(tpm)];
        if tpm <= 3 && dartel
            job.tissue(tpm).native=[1 1];
        end
        if tpm >= 4
            job.tissue(tpm).native=[0 0];
        end
    end
    
    job.tissue(tpm+1:end)=[];
    if del
        job.warp.write=[0,0];
    else
        job.warp.write=[1,1];
    end
    
    spm_preproc_run(job); % run "Segment" in SPM 12 (Old "Segment" is now referred to as "Old Segment").
    % delete unused files
    [~,fn]=fileparts(files{1});
    ea_delete([directory,fn,'_seg8.mat']);
    
    disp('Done.');
end





