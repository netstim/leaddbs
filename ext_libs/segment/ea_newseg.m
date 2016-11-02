function ea_newseg(directory,file,dartel,options,del)
% SPM NewSegment
% we cannot generate the TPM from the SPM TPM anymore since
% we use the enhanced TPM by Lorio / Draganski:
% http://unil.ch/lren/home/menuinst/data--utilities.html

if nargin < 5
    del = 1;
end

if ~dartel && exist([directory, 'c1', file], 'file') || dartel && exist([directory, 'rc1', file], 'file')
    disp('Segmentation already done!');
else
    disp('Segmentation...');

    load([options.earoot,'ext_libs',filesep,'segment',filesep,'segjob12']);
    tpminf=[options.earoot,'templates',filesep,'TPM_Lorio_Draganski.nii'];
    job.channel.vols{1}=[directory,file,',1'];
    for tpm=1:6
        job.tissue(tpm).tpm=[tpminf,',',num2str(tpm)];
        if tpm <= 3 && dartel
            job.tissue(tpm).native=[1 1];
        end
        if tpm >= 4
            job.tissue(tpm).native=[0 0];
        end
    end
    job.warp.write=[1,1];
    spm_preproc_run(job); % run "Segment" in SPM 12 (Old "Segment" is now referred to as "Old Segment").

    % delete unused files
    [~,fn]=fileparts(file);
    ea_delete([directory,fn,'_seg8.mat']);

    % del==0: keep the deformation field when for ea_normalize_spmnewseg
    if del
        ea_delete([directory, 'y_', file]);
        ea_delete([directory, 'iy_', file]);
    end

    disp('Done.');
end
