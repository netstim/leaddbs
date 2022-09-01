function ea_create_tpm_darteltemplate(mute)

if exist([ea_space,'TPM.nii'],'file') && ~ea_checktpmresolution && exist([ea_space,'.tpm'],'file')
    return
end
if ~exist('mute','var')
    answ=questdlg('Lead Neuroimaging Suite needs to generate some files needed for the process. This only needs to be done once but will take some additional time. The process you started will be performed afterwards.','Additional files needed','Proceed','Abort','Proceed');
    switch answ
        case 'Abort'
            ea_error('Process aborted by user');
    end
end

load([ea_space,'spacedef.mat']);

if ~exist([ea_space,'dartel'], 'dir')
    mkdir([ea_space,'dartel']);
end

% TPM Lorio Draganski to be kept in 'templates' folder since will be used to generate TPM in each space.
tpmfile=[ea_getearoot,'templates',filesep,'TPM_Lorio_Draganski.nii'];

rebuildtpm=1;
if isfield(spacedef,'tpm')

    switch spacedef.tpm
        case 'custom_fixed'
            tpmfile=[ea_space,'TPM.nii'];
            rebuildtpm=0;
        case 'custom_rebuild'
            tpmfile=[ea_space,'TPM.nii'];
            rebuildtpm=1;

            ea_backuprestore([ea_space,'TPM.nii']);
    end
end

tpmHdr = ea_open_vol(tpmfile);
tpmnum = tpmHdr.volnum;

[~,ix]=ismember({'STN','GPi','GPe','RN'},spacedef.templates);
if any(ix)
    spacedef.templates(ix)=[];
end

if rebuildtpm
    for t=1:length(spacedef.templates)
        matlabbatch{1}.spm.spatial.preproc.channel(t).vols = {[ea_space,spacedef.templates{t},'.nii,1']};
        matlabbatch{1}.spm.spatial.preproc.channel(t).biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel(t).biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel(t).write = [0 0];
    end

    for tpspectr=1:tpmnum
        matlabbatch{1}.spm.spatial.preproc.tissue(tpspectr).tpm = {[tpmfile,',',num2str(tpspectr)]};
        matlabbatch{1}.spm.spatial.preproc.tissue(tpspectr).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(tpspectr).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(tpspectr).warped = [0 0];
    end

    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];

    spm_jobman('run',{matlabbatch});
    clear matlabbatch

    ea_delete([ea_space,spacedef.templates{1},'_seg8.mat']);
else
    % split TPM
    matlabbatch{1}.spm.util.split.vol = {[ea_space,'TPM.nii,1']};
    matlabbatch{1}.spm.util.split.outdir = {ea_space};
    spm_jobman('run',{matlabbatch});
    clear matlabbatch

    for t=1:tpmnum
       movefile([ea_space,'TPM_',sprintf('%05.f',t),'.nii'],[ea_space,'c',num2str(t),spacedef.templates{1},'.nii']);
    end
end

for c=1:tpmnum
    if c<3
       copyfile([ea_space,'c',num2str(c),spacedef.templates{1},'.nii'],[ea_space,'c',num2str(c),'mask.nii']);
    end
    movefile([ea_space,'c',num2str(c),spacedef.templates{1},'.nii'],[ea_space([],'dartel'),filesep,'dartelmni_6_hires_',sprintf('%05d',c),'.nii']);
end

prefs=ea_prefs('');
if ~startsWith(spacedef.tpm,'custom')
    for c=1:tpmnum
        fina=[ea_space([],'dartel'),'dartelmni_6_hires_',sprintf('%05d',c),'.nii'];
        nii=ea_load_nii(fina); % change datatype to something high for reslicing and smoothing.
        nii.dt(1) = 16;
        ea_write_nii(nii);
        if ~(prefs.normalize.spm.resolution==0.5) % reslice images
            ea_reslice_nii(fina,fina,ones(1,3)*prefs.normalize.spm.resolution,1,[],1);
        end
        % apply very light smooth on custom TPMs
        job{1}.spm.spatial.smooth.data = {fina};
        job{1}.spm.spatial.smooth.fwhm = [0.5 0.5 0.5];
        job{1}.spm.spatial.smooth.dtype = 0;
        job{1}.spm.spatial.smooth.im = 0;
        job{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',{job});
        clear job
        [pth,fn,ext]=fileparts(fina);
        movefile(fullfile(pth,['s',fn,ext]),fullfile(pth,[fn,ext]));

        nii=ea_load_nii(fina); % change datatype back to uint8
        nii.dt(1) = 2;
        ea_delete(fina);
        ea_write_nii(nii);

        matlabbatch{1}.spm.util.cat.vols{c} = fina;
    end

    matlabbatch{1}.spm.util.cat.vols = matlabbatch{1}.spm.util.cat.vols';

    matlabbatch{1}.spm.util.cat.name = [ea_space,'TPM.nii'];
    matlabbatch{1}.spm.util.cat.dtype = 16;
    spm_jobman('run',{matlabbatch});
    clear matlabbatch
    ea_delete([ea_space,'TPM.mat']);
end

% make sure TPM sums to 1 everywhere
nii=ea_load_untouch_nii([ea_space,'TPM.nii']);
nii.img=nii.img./repmat(sum(nii.img,4),1,1,1,tpmnum);
ea_save_untouch_nii(nii,[ea_space,'TPM.nii']);

wd=ea_space([],'dartel');
% gunzip([wd,'dartelmni_6_hires.nii.gz']);
% spm_file_split([wd,'dartelmni_6_hires.nii']);
gs=[0,2,3,5,6,8];
expo=6:-1:1;
spacedef=ea_getspacedef;
if isfield(spacedef,'dartelpreserve')
    dpres=spacedef.dartelpreserve;
else
    dpres=3;
end

for s=1:6
    for tpm=1:dpres
        % smooth
        if gs(s)
            matlabbatch{1}.spm.spatial.smooth.data = {[wd,'dartelmni_6_hires_',sprintf('%05d',tpm),'.nii,1']};
            matlabbatch{1}.spm.spatial.smooth.fwhm = [gs(s),gs(s),gs(s)];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = ['s',num2str(gs(s))];
            jobs{1}=matlabbatch;
            spm_jobman('run',jobs);
            clear jobs matlabbatch
        else
            copyfile([wd,'dartelmni_6_hires_',sprintf('%05d',tpm),'.nii'],[wd,'s0','dartelmni_6_hires_',sprintf('%05d',tpm),'.nii']);
        end
        clear jobs matlabbatch

        % set to resolution of TPM file
        matlabbatch{1}.spm.util.imcalc.input = {[ea_space,'TPM.nii,1'];
            [wd,'s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',tpm),'.nii']};
        matlabbatch{1}.spm.util.imcalc.output = [wd,'s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',tpm),'.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {wd};
        matlabbatch{1}.spm.util.imcalc.expression = 'i2';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 5;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
        jobs{1}=matlabbatch;
        spm_jobman('run',jobs);
        clear jobs matlabbatch
    end

    for tpm=1:dpres
        matlabbatch{1}.spm.util.cat.vols{tpm} = [wd,'s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',tpm),'.nii'];
    end
    matlabbatch{1}.spm.util.cat.vols=matlabbatch{1}.spm.util.cat.vols';
    matlabbatch{1}.spm.util.cat.name = [wd,'dartelmni_',num2str(expo(s)),'.nii'];
    matlabbatch{1}.spm.util.cat.dtype = 0;
    jobs{1}=matlabbatch;
    spm_jobman('run',jobs);
    clear jobs matlabbatch

    disp('Cleaning up.');

    % cleanup
    ea_delete([wd,'s',num2str(gs(s)),'dartelmni_6_hires_00*.*']);
end

% further cleanup
ea_delete([wd,'dartelmni_*.mat']);
for c=1:tpmnum
    ea_delete([wd,'dartelmni_6_hires_',sprintf('%05d',c),'.nii']);
end
% gzip([wd,'dartelmni_6_hires.nii']);
% ea_delete([wd,'dartelmni_6_hires.nii']);
disp('Done.');

ea_addshoot;

% %Finalize TPM
% % Split TPM:
%         matlabbatch{1}.spm.util.split.vol = {[ea_space,'TPM.nii,1']};
%         matlabbatch{1}.spm.util.split.outdir = {ea_space};
%         spm_jobman('run',{matlabbatch});
%         clear matlabbatch
%
%         copyfile([ea_space,'atlas.nii'],[ea_space,'catlas.nii']);
%         ea_conformspaceto([ea_space,'TPM_',sprintf('%05d',1),'.nii'],[ea_space,'catlas.nii'],1);
%
%         c1=ea_load_nii([ea_space,'TPM_',sprintf('%05d',1),'.nii']);
%         if exist([ea_space,'atlas.nii'],'file') % add atlas.
%
%             atlas=ea_load_nii([ea_space,'catlas.nii']);
%             c1.img(atlas.img>0.1)=atlas.img(atlas.img>0.1);
%             ea_write_nii(c1);
%             c2=ea_load_nii([ea_space,'TPM_',sprintf('%05d',2),'.nii']);
%             c2.img(atlas.img>0.1)=0;
%             ea_write_nii(c2);
%             c3=ea_load_nii([ea_space,'TPM_',sprintf('%05d',3),'.nii']);
%             c3.img(atlas.img>0.1)=0;
%             ea_write_nii(c3);
%
%         end
%
%         matlabbatch{1}.spm.spatial.smooth.data = {
%             [ea_space,'TPM_00001.nii,1']
%             [ea_space,'TPM_00002.nii,1']
%             [ea_space,'TPM_00003.nii,1']
%             [ea_space,'TPM_00004.nii,1']
%             [ea_space,'TPM_00005.nii,1']
%             [ea_space,'TPM_00006.nii,1']
%             };
%         matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
%         matlabbatch{1}.spm.spatial.smooth.dtype = 0;
%         matlabbatch{1}.spm.spatial.smooth.im = 0;
%         matlabbatch{1}.spm.spatial.smooth.prefix = 's';
%         spm_jobman('run',{matlabbatch}); clear matlabbatch
%
%         matlabbatch{1}.spm.util.cat.vols = {
%             [ea_space,'sTPM_00001.nii,1']
%             [ea_space,'sTPM_00002.nii,1']
%             [ea_space,'sTPM_00003.nii,1']
%             [ea_space,'sTPM_00004.nii,1']
%             [ea_space,'sTPM_00005.nii,1']
%             [ea_space,'sTPM_00006.nii,1']
%             };
%         matlabbatch{1}.spm.util.cat.name = [ea_space,'TPM.nii'];
%         matlabbatch{1}.spm.util.cat.dtype = 16;
%     spm_jobman('run',{matlabbatch}); clear matlabbatch
% % make sure TPM sums to 1 everywhere
%
% nii=ea_load_untouch_nii([ea_space,'TPM.nii']);
% nii.img=nii.img./repmat(sum(nii.img,4),1,1,1,tpmnum);
% ea_save_untouch_nii(nii,[ea_space,'TPM.nii']);
%
% ea_delete([ea_space,'*PM_0*.nii']);
% ea_delete([ea_space,'TPM.mat']);
% try
%     ea_delete([ea_space,'catlas.nii']);
% end

% set marker that last gen TPM has been set.
fid=fopen([ea_space,'.tpm'],'w');
fprintf(fid,'%s','Novel TPM generated by Lead Neuroimaging Suite.');
fclose(fid);


function ea_addshoot
spacedef=ea_getspacedef;
if isfield(spacedef,'dartelpreserve')
    dpres=spacedef.dartelpreserve;
else
    dpres=3;
end

%if ~exist([ea_space([],'dartel'),filesep,'shootmni_1.nii'],'file');
    root=[ea_space([],'dartel'),filesep];
    for dt=1:6
        nii=ea_load_nii([root,'dartelmni_',num2str(dt),'.nii']);

        matlabbatch{1}.spm.util.split.vol = {[root,'dartelmni_',num2str(dt),'.nii,1']};
        matlabbatch{1}.spm.util.split.outdir = {root};
        spm_jobman('run',{matlabbatch});
        clear matlabbatch

        X=-sum(nii.img,4);
        X=X+1;

        nii=ea_load_nii([root,'dartelmni_',num2str(dt),'_',sprintf('%05.0f',1),'.nii']);
        nii.img=X;
        nii.fname=[root,'/dartelmni_',num2str(dt),'_',sprintf('%05.0f',dpres+1),'.nii'];
        ea_write_nii(nii);

        for i=1:dpres+1
            matlabbatch{1}.spm.util.cat.vols{i} = [root,'dartelmni_',num2str(dt),'_',sprintf('%05.0f',i),'.nii'];
        end
        matlabbatch{1}.spm.util.cat.vols=matlabbatch{1}.spm.util.cat.vols';
        matlabbatch{1}.spm.util.cat.name = ['shootmni_',num2str(dt),'.nii'];
        matlabbatch{1}.spm.util.cat.dtype = 0;
        spm_jobman('run',{matlabbatch});
        clear matlabbatch

        for i=1:dpres+1 % cleanup
            ea_delete([root,'dartelmni_',num2str(dt),'_',sprintf('%05.0f',i),'.nii']);
        end
        ea_delete([root,'shootmni_',num2str(dt),'.mat']);
    end
%end
