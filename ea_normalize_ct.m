function varargout=ea_normalize_ct(options)



if ischar(options) % return name of method.
    varargout{1}='Fuse CT and MRI [Not robust]';
    return
end


if ~exist([options.earoot,'templates',filesep,'TPM.nii'],'file')
   ea_generate_tpm; 
    
end

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



%% apply segmentations

matlabbatch{1}.spm.util.defs.comp{1}.def = {[root,pt,filesep,'y_',moving]};
matlabbatch{1}.spm.util.defs.ofname = '';
matlabbatch{1}.spm.util.defs.fnames = {[root,pt,filesep,moving]};
matlabbatch{1}.spm.util.defs.savedir.savedef = 1;
matlabbatch{1}.spm.util.defs.interp = 1;

jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs

matlabbatch{1}.spm.util.defs.comp{1}.def = {[root,pt,filesep,'y_',fixed]};
matlabbatch{1}.spm.util.defs.ofname = '';
matlabbatch{1}.spm.util.defs.fnames = {[root,pt,filesep,fixed]};
matlabbatch{1}.spm.util.defs.savedir.savedef = 1;
matlabbatch{1}.spm.util.defs.interp = 1;

jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs



%% Coreg C1 volumes Pre to CT
try

matlabbatch{1}.spm.spatial.normalise.est.subj.source = {[options.root,options.prefs.patientdir,filesep,'w',options.prefs.rawctnii_unnormalized,',1']};
matlabbatch{1}.spm.spatial.normalise.est.subj.wtsrc = {[options.root,options.prefs.patientdir,filesep,'wc1',options.prefs.rawctnii_unnormalized,',1']};
matlabbatch{1}.spm.spatial.normalise.est.eoptions.template = {[options.root,options.prefs.patientdir,filesep,'w',options.prefs.prenii_unnormalized,',1']};
matlabbatch{1}.spm.spatial.normalise.est.eoptions.weight = {[options.root,options.prefs.patientdir,filesep,'wc1',options.prefs.prenii_unnormalized,',1']};
matlabbatch{1}.spm.spatial.normalise.est.eoptions.smosrc = 10;
matlabbatch{1}.spm.spatial.normalise.est.eoptions.smoref = 10;
matlabbatch{1}.spm.spatial.normalise.est.eoptions.regtype = 'mni';
matlabbatch{1}.spm.spatial.normalise.est.eoptions.cutoff = 25;
matlabbatch{1}.spm.spatial.normalise.est.eoptions.nits = 0;
matlabbatch{1}.spm.spatial.normalise.est.eoptions.reg = 1;

jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs




% apply estimated transformations to CT image.

    
    [~,nm]=fileparts(['w',options.prefs.rawctnii_unnormalized]); % cut off file extension

    voxi=[1 1 1]; % export highres
    bbi=nan(2,3);
    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname = {[options.root,options.prefs.patientdir,filesep,nm,'_sn.mat']};
    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.vox = voxi;
    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.bb = bbi;
    
    matlabbatch{1}.spm.util.defs.ofname = 'ea_normparams';
    matlabbatch{1}.spm.util.defs.fnames = {[options.root,options.prefs.patientdir,filesep,'w',options.prefs.rawctnii_unnormalized,',1']}; % fusion.nii
    matlabbatch{1}.spm.util.defs.savedir.saveusr = {[options.root,options.prefs.patientdir,filesep]};
    matlabbatch{1}.spm.util.defs.interp = 6;
    
    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear matlabbatch jobs;


%% call up checkreg

matlabbatch{1}.spm.util.checkreg.data = {
    [options.root,options.prefs.patientdir,filesep,'ww',options.prefs.rawctnii_unnormalized,',1']
    [options.root,options.prefs.patientdir,filesep,'w',options.prefs.prenii_unnormalized,',1']
    };

jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs

%% Fuse images:

    
    matlabbatch{1}.spm.util.imcalc.input = {[options.root,options.prefs.patientdir,filesep,'w',options.prefs.prenii_unnormalized,',1']
        [options.root,options.prefs.patientdir,filesep,'ww',options.prefs.rawctnii_unnormalized,',1']};
    matlabbatch{1}.spm.util.imcalc.output = ['ww',options.prefs.ctnii_unnormalized]; % usually fusion.
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


%% Fuse images:

    
    matlabbatch{1}.spm.util.imcalc.input = {[options.root,options.prefs.patientdir,filesep,'w',options.prefs.prenii_unnormalized,',1']
        [options.root,options.prefs.patientdir,filesep,'w',options.prefs.rawctnii_unnormalized,',1']};
    matlabbatch{1}.spm.util.imcalc.output = [options.prefs.ctnii_unnormalized]; % usually fusion.
    matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.prefs.patientdir,filesep]};
    matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2)/2';
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = -2;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    
    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear matlabbatch jobs















% make normalization "permanent" and include correct bounding box.


outf=options.prefs.ctnii;
fina=[options.root,options.prefs.patientdir,filesep,options.prefs.ctnii_unnormalized,',1'];


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



% rename to lfusion.nii
movefile([options.root,options.prefs.patientdir,filesep,options.prefs.ctnii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.gctnii]);




% cleanup

delete([options.root,options.prefs.patientdir,filesep,'wc1',options.prefs.prenii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'wc2',options.prefs.prenii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'wc3',options.prefs.prenii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'wc4',options.prefs.prenii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'c1',options.prefs.prenii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'c2',options.prefs.prenii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'c3',options.prefs.prenii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'c4',options.prefs.prenii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'c5',options.prefs.prenii_unnormalized]);

delete([options.root,options.prefs.patientdir,filesep,'wc1',options.prefs.rawctnii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'wc2',options.prefs.rawctnii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'wc3',options.prefs.rawctnii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'wc4',options.prefs.rawctnii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'c1',options.prefs.rawctnii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'c2',options.prefs.rawctnii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'c3',options.prefs.rawctnii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'c4',options.prefs.rawctnii_unnormalized]);
delete([options.root,options.prefs.patientdir,filesep,'c5',options.prefs.rawctnii_unnormalized]);


function ea_segment(path)

spmdir=spm('dir');
matlabbatch{1}.spm.tools.preproc8.channel.vols = {[path,',1']};
matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[spmdir,filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,1']};
matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[spmdir,filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,2']};
matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[spmdir,filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,3']};
matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[spmdir,filesep,'toolbox',filesep,'Seg',filesep,'TPM.nii,4']};
matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [1 0];
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
matlabbatch{1}.spm.tools.preproc8.warp.write = [0 1];

jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs