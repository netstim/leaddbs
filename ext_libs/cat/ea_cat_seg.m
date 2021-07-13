function ea_cat_seg(options)


options.prefs=ea_prefs;
[options,presentfiles]=ea_assignpretra(options);

job.data = {[options.root,options.patientname,filesep,presentfiles{1}]};
job.nproc = 0;
job.opts.tpm = {[ea_space,'TPM.nii']};
job.opts.affreg = 'mni';
job.opts.biasstr = 0.5;
job.extopts.APP = 1070;
job.extopts.LASstr = 0.5;
job.extopts.gcutstr = 0.5;
job.extopts.cleanupstr = 0.5;
job.extopts.registration.darteltpm = {[ea_space,'dartel',filesep,'dartelmni_1.nii']};
for tp=1:6
    job.extopts.darteltpms{tp}=[ea_space,'dartel',filesep,'dartelmni_',num2str(tp),'.nii'];
end
job.extopts.registration.shootingtpm = {[ea_space,'dartel',filesep,'shootmni_1.nii']};
for tp=1:6
    job.extopts.shootingtpms{tp}=[ea_space,'dartel',filesep,'shootmni_',num2str(tp),'.nii'];
end
job.extopts.registration.regstr = 0.5;
job.extopts.vox = 1.5;
job.output.surface = 1;
job.output.ROI = 0;
job.output.GM.native = 0;
job.output.GM.mod = 0;
job.output.GM.dartel = 0;
job.output.WM.native = 0;
job.output.WM.mod = 0;
job.output.WM.dartel = 0;
job.output.bias.warped = 0;
job.output.jacobian.warped = 0;
job.output.warps = [0 0];

% run CAT12 surface reconstruction:
ea_cat_run(job);

% visualize results using surfice:
script=['BEGIN; ',...
    'RESETDEFAULTS; ',...
    'MESHLOAD(''',ea_path_helper([options.root,options.patientname,filesep,'surf',filesep,'lh.central.',ea_stripext(presentfiles{1}),'.gii']),'''); ',...
    'OVERLAYLOAD(''',ea_path_helper([options.root,options.patientname,filesep,'surf',filesep,'lh.thickness.',ea_stripext(presentfiles{1})]),'''); ',...
    ' OVERLAYCOLORNAME(1, ''ACTC'');',...
    ' OVERLAYMINMAX(1,0,4);',...
    'END.'];

ea_surfice(script,0);
