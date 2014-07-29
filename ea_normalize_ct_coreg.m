function ea_normalize_ct_coreg(options)


alsouseregutil=0;

if exist([options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,'.gz'],'file')
    try
        gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,'.gz']);
    end
    try
        gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.rawctnii_unnormalized,'.gz']);
    end
    
end



% set local file variables
root=options.root;
pt=options.prefs.patientdir;

moving=options.prefs.rawctnii_unnormalized;
fixed=options.prefs.prenii_unnormalized;


%% segment both CT and Preoperative MR using SPM's New Segment routine
% (Segment in SPM12).
ea_segment([root,pt,filesep,moving]);
ea_segment([root,pt,filesep,fixed]);




%% Coreg C1 volumes Pre to CT
copyfile([options.root,options.prefs.patientdir,filesep,'c1',options.prefs.rawctnii_unnormalized],[options.root,options.prefs.patientdir,filesep,'kc1',options.prefs.rawctnii_unnormalized]);
copyfile([options.root,options.prefs.patientdir,filesep,options.prefs.rawctnii_unnormalized],[options.root,options.prefs.patientdir,filesep,'k',options.prefs.rawctnii_unnormalized]);

matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[options.root,options.prefs.patientdir,filesep,'c1',options.prefs.prenii_unnormalized,',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[options.root,options.prefs.patientdir,filesep,'kc1',options.prefs.rawctnii_unnormalized,',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.other = {[options.root,options.prefs.patientdir,filesep,'k',options.prefs.rawctnii_unnormalized,',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];


jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs


%% call up checkreg

matlabbatch{1}.spm.util.checkreg.data = {
    [options.root,options.prefs.patientdir,filesep,'k',options.prefs.rawctnii_unnormalized,',1']
    [options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1']
    };

jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs




if alsouseregutil

    
movnii=load_untouch_nii([options.root,options.prefs.patientdir,filesep,options.prefs.rawctnii_unnormalized]);
mask_movnii=load_untouch_nii([options.root,options.prefs.patientdir,filesep,'c1',options.prefs.rawctnii_unnormalized]);
fixnii=load_untouch_nii([options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized]);
mask_fixnii=load_untouch_nii([options.root,options.prefs.patientdir,filesep,'c1',options.prefs.prenii_unnormalized]);

regopts.Registration='affine';
regopts.Similarity='mi';
%regopts.MaskMoving=double(mask_movnii.img);
%regopts.MaskStatic=double(mask_fixnii.img);

[Ireg] = image_registration(double(movnii.img),double(fixnii.img),regopts);
fixnii.img=Ireg;
save_untouch_nii(fixnii,[options.root,options.prefs.patientdir,filesep,'rk',options.prefs.rawctnii_unnormalized])

matlabbatch{1}.spm.util.checkreg.data = {
    [options.root,options.prefs.patientdir,filesep,'rk',options.prefs.rawctnii_unnormalized,',1']
    [options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1']
    };

jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs
end



%% Fuse images:
for fusion=1:(1+alsouseregutil) % repeat if also use regutil is set on.
    
    if fusion==2
        addstr='r';
    else
        addstr='';
    end
    
    matlabbatch{1}.spm.util.imcalc.input = {[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1']
        [options.root,options.prefs.patientdir,filesep,addstr,'k',options.prefs.rawctnii_unnormalized,',1']};
    matlabbatch{1}.spm.util.imcalc.output = [addstr,options.prefs.ctnii_unnormalized]; % usually fusion.
    matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.prefs.patientdir,filesep]};
    matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2)/2';
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = -2;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    
    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear matlabbatch jobs
end






%% segment preop MRI

    matlabbatch{1}.spm.spatial.preproc.data = {[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1']};
    matlabbatch{1}.spm.spatial.preproc.output.GM = [0 0 1];
    matlabbatch{1}.spm.spatial.preproc.output.WM = [0 0 1];
    matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 0 0];
    matlabbatch{1}.spm.spatial.preproc.output.biascor = 0;
    matlabbatch{1}.spm.spatial.preproc.output.cleanup = 0;
    matlabbatch{1}.spm.spatial.preproc.opts.tpm = {
        fullfile(options.earoot,'templates','mni_icbm152_gm_tal_nlin_asym_09c.nii')
        fullfile(options.earoot,'templates','mni_icbm152_wm_tal_nlin_asym_09c.nii')
        fullfile(options.earoot,'templates','mni_icbm152_csf_tal_nlin_asym_09c.nii')
        };
    matlabbatch{1}.spm.spatial.preproc.opts.ngaus = [2
        2
        2
        4];
    matlabbatch{1}.spm.spatial.preproc.opts.regtype = 'mni'; %'mni';
    matlabbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
    matlabbatch{1}.spm.spatial.preproc.opts.warpco = 25;
    matlabbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.opts.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.opts.msk = {''};

jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs




% apply estimated transformations to fusion image.
for fusion=1:(1+alsouseregutil) % repeat if also use regutil is set on.
    
    if fusion==2
        addstr='r';
    else
        addstr='';
    end
    
    [~,nm]=fileparts(options.prefs.prenii_unnormalized); % cut off file extension
    
    voxi=[0.22 0.22 0.5]; % export highres
    bbi=[-55 45 9.5; 55 -65 -25]; % with small bounding box
    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname = {[options.root,options.prefs.patientdir,filesep,nm,'_seg_sn.mat']};
    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.vox = voxi;
    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.bb = bbi;
    
    matlabbatch{1}.spm.util.defs.ofname = 'ea_normparams';
    matlabbatch{1}.spm.util.defs.fnames = {[options.root,options.prefs.patientdir,filesep,addstr,options.prefs.ctnii_unnormalized,',1']}; % fusion.nii
    matlabbatch{1}.spm.util.defs.savedir.saveusr = {[options.root,options.prefs.patientdir,filesep]};
    matlabbatch{1}.spm.util.defs.interp = 6;
    
    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear matlabbatch jobs;
    
end


% make normalization "permanent" and include correct bounding box.


outf=options.prefs.ctnii;
fina=[options.root,options.prefs.patientdir,filesep,'w',options.prefs.ctnii_unnormalized,',1'];


% save a backup

if exist([options.root,options.prefs.patientdir,filesep,outf],'file')
    movefile([options.root,options.prefs.patientdir,filesep,outf],[options.root,options.prefs.patientdir,filesep,'ea_backup',date,num2str(now),outf]);
end

matlabbatch{1}.spm.util.imcalc.input = {[options.earoot,'templates',filesep,'bb.nii,1']
    fina
    };
matlabbatch{1}.spm.util.imcalc.output = outf;
matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.prefs.patientdir,filesep]};
matlabbatch{1}.spm.util.imcalc.expression = ['i2'];
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
jobs{1}=matlabbatch;
    cfg_util('run',jobs);

clear matlabbatch jobs;




% build global versions of files



[~,nm]=fileparts(options.prefs.prenii_unnormalized); % cut off file extension

voxi=[0.5 0.5 0.5]; % export highres
bbi=nan(2,3);
matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname = {[options.root,options.prefs.patientdir,filesep,nm,'_seg_sn.mat']};
matlabbatch{1}.spm.util.defs.comp{1}.sn2def.vox = voxi;
matlabbatch{1}.spm.util.defs.comp{1}.sn2def.bb = bbi;

matlabbatch{1}.spm.util.defs.ofname = 'ea_normparams';
matlabbatch{1}.spm.util.defs.fnames = {[options.root,options.prefs.patientdir,filesep,options.prefs.ctnii_unnormalized,',1']}; % fusion.nii
matlabbatch{1}.spm.util.defs.savedir.saveusr = {[options.root,options.prefs.patientdir,filesep]};
matlabbatch{1}.spm.util.defs.interp = 6;

jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs;

% rename to lfusion.nii
movefile([options.root,options.prefs.patientdir,filesep,'w',options.prefs.ctnii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.gctnii]);




function ea_segment(path)

spmdir=spm('dir');
matlabbatch{1}.spm.tools.preproc8.channel.vols = {[path,',1']};
matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[spmdir,filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,1']};
matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[spmdir,filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,2']};
matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[spmdir,filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,3']};
matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[spmdir,filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,4']};
matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[spmdir,filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,5']};
matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[spmdir,filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,6']};
matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.warp.mrf = 0;
matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
matlabbatch{1}.spm.tools.preproc8.warp.write = [0 0];

jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs