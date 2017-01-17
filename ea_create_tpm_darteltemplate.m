function ea_create_tpm_darteltemplate(mute)

if exist([ea_space,'TPM.nii'],'file') && ~ea_checktpmresolution && exist([ea_space,'.newtpm'],'file')
    return
end
if ~exist('mute','var')
    answ=questdlg('Lead Neuroimaging Suite needs to generate some files needed for the process. This only needs to be done once but will take some additional time. The process you started will be performed afterwards.','Additional files needed','Proceed','Abort','Proceed');
    switch answ
        case 'Abort'
            ea_error('Process aborted by user');
    end
end

load([ea_space,'ea_space_def.mat']);
for t=1:length(spacedef.templates) 
    matlabbatch{1}.spm.spatial.preproc.channel(t).vols = {[ea_space,spacedef.templates{t},'.nii,1']};
    matlabbatch{1}.spm.spatial.preproc.channel(t).biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel(t).biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel(t).write = [0 0];
end
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[ea_getearoot,'templates',filesep,'TPM_Lorio_Draganski.nii,1']}; % This is correct ? TPM Lorio Draganski to be kept in /templates folder since will be used to generate TPM in each space.
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[ea_getearoot,'templates',filesep,'TPM_Lorio_Draganski.nii,2']};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[ea_getearoot,'templates',filesep,'TPM_Lorio_Draganski.nii,3']};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[ea_getearoot,'templates',filesep,'TPM_Lorio_Draganski.nii,4']};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[ea_getearoot,'templates',filesep,'TPM_Lorio_Draganski.nii,5']};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[ea_getearoot,'templates',filesep,'TPM_Lorio_Draganski.nii,6']};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];

spm_jobman('run',{matlabbatch});
clear matlabbatch

delete([ea_space,spacedef.templates{1},'_seg8.mat']); 
if ~exist([ea_space,'dartel'], 'dir')
    mkdir([ea_space,'dartel']);
end
for c=1:6
    if c<3
       copyfile([ea_space,'c',num2str(c),spacedef.templates{1},'.nii'],[ea_space,'c',num2str(c),'mask.nii']);
    end
    movefile([ea_space,'c',num2str(c),spacedef.templates{1},'.nii'],[ea_space([],'dartel'),filesep,'dartelmni_6_hires_',sprintf('%05d',c),'.nii']);
end

% add atlas
if exist([ea_space,'atlas.nii'],'file')
    copyfile([ea_space,'atlas.nii'],[ea_space([],'dartel'),filesep,'atlas.nii']);
    ea_conformspaceto([ea_space([],'dartel'),filesep,'dartelmni_6_hires_',sprintf('%05d',1),'.nii'],[ea_space([],'dartel'),filesep,'atlas.nii'],6);
    
    c1=ea_load_nii([ea_space([],'dartel'),filesep,'dartelmni_6_hires_',sprintf('%05d',1),'.nii']);
    atlas=ea_load_nii([ea_space([],'dartel'),filesep,'atlas.nii']);
    c1.img(atlas.img>0.1)=atlas.img(atlas.img>0.1);
    ea_write_nii(c1);
    c2=ea_load_nii([ea_space([],'dartel'),filesep,'dartelmni_6_hires_',sprintf('%05d',2),'.nii']);
    c2.img(atlas.img>0.1)=0;
    ea_write_nii(c2);
    c3=ea_load_nii([ea_space([],'dartel'),filesep,'dartelmni_6_hires_',sprintf('%05d',3),'.nii']);
    c3.img(atlas.img>0.1)=0;
    ea_write_nii(c3);
end
prefs=ea_prefs('');

for c=1:6
    fina=[ea_space([],'dartel'),filesep,'dartelmni_6_hires_',sprintf('%05d',c),'.nii'];
    nii=ea_load_nii(fina); % change datatype to something high for reslicing and smoothing.
    nii.dt=[16,0];
    ea_write_nii(nii);
    if ~(prefs.normalize.spm.resolution==0.5) % reslice images
        ea_reslice_nii(fina,fina,[prefs.normalize.spm.resolution prefs.normalize.spm.resolution prefs.normalize.spm.resolution],1,[],6);
    end
    % apply very light smooth
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
    nii.dt=[2,0];
    delete(fina);
    ea_write_nii(nii);
    matlabbatch{1}.spm.util.cat.vols{c} = fina;
end

matlabbatch{1}.spm.util.cat.vols = matlabbatch{1}.spm.util.cat.vols';

matlabbatch{1}.spm.util.cat.name = [ea_space,'TPM.nii'];
matlabbatch{1}.spm.util.cat.dtype = 16;
spm_jobman('run',{matlabbatch});
clear matlabbatch

delete([ea_space,'TPM.mat']);

wd=[ea_space([],'dartel'),filesep];
%gunzip([wd,'dartelmni_6_hires.nii.gz']);
%spm_file_split([wd,'dartelmni_6_hires.nii']);
gs=[0,2,3,5,6,8];
expo=6:-1:1;
for s=1:6
    for tpm=1:3 
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
    
    matlabbatch{1}.spm.util.cat.vols = {[wd,'s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',1),'.nii'];
        [wd,'s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',2),'.nii'];
        [wd,'s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',3),'.nii']};
    matlabbatch{1}.spm.util.cat.name = [wd,'dartelmni_',num2str(expo(s)),'.nii'];
    matlabbatch{1}.spm.util.cat.dtype = 0;
    jobs{1}=matlabbatch;
    spm_jobman('run',jobs);
    clear jobs matlabbatch
    
    disp('Cleaning up.');
    
    % cleanup
    delete([wd,'s',num2str(gs(s)),'dartelmni_6_hires_00*.*']);
end

% further cleanup
delete([wd,'dartelmni_*.mat']);
for c=1:6
    delete([wd,'dartelmni_6_hires_',sprintf('%05d',c),'.nii']);
end
%gzip([wd,'dartelmni_6_hires.nii']);
%delete([wd,'dartelmni_6_hires.nii']);
disp('Done.');

ea_addshoot;


% lightly smooth TPM
for tp=1:6
matlabbatch{1}.spm.spatial.smooth.data{tp} = [ea_space,'TPM.nii,',num2str(tp)];
end
matlabbatch{1}.spm.spatial.smooth.data=matlabbatch{1}.spm.spatial.smooth.data';
matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',{matlabbatch});
clear matlabbatch
movefile([ea_space,'sTPM.nii'],[ea_space,'TPM.nii']);
% set marker that last gen TPM has been set.
fid=fopen([ea_space,'.newtpm'],'w');
fprintf(fid,'%s','Novel TPM generated by Lead Neuroimaging Suite.');
fclose(fid);


function ea_addshoot
if ~exist([ea_space([],'dartel'),filesep,'shootmni_1.nii'],'file');
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
        nii.fname=[root,'/dartelmni_',num2str(dt),'_',sprintf('%05.0f',4),'.nii'];
        ea_write_nii(nii);

        for i=1:4
            matlabbatch{1}.spm.util.cat.vols{i} = [root,'dartelmni_',num2str(dt),'_',sprintf('%05.0f',i),'.nii'];
        end
        matlabbatch{1}.spm.util.cat.vols=matlabbatch{1}.spm.util.cat.vols';
        matlabbatch{1}.spm.util.cat.name = ['shootmni_',num2str(dt),'.nii'];
        matlabbatch{1}.spm.util.cat.dtype = 0;
        spm_jobman('run',{matlabbatch});
        clear matlabbatch

        for i=1:4 % cleanup
            delete([root,'dartelmni_',num2str(dt),'_',sprintf('%05.0f',i),'.nii']);
        end
        delete([root,'shootmni_',num2str(dt),'.mat']);
    end
end
