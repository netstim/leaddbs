function    ea_checkfiles(options)
% load files
if strcmp(options.prefs.patientdir,'Choose Patient Directory')
    ea_error('Please choose patient directory first');
end


tranii=load_untouch_nii([options.root,options.prefs.patientdir,filesep,options.prefs.tranii]);
try
cornii=load_untouch_nii([options.root,options.prefs.patientdir,filesep,options.prefs.cornii]);
end
% check voxel dimensions of each file
exp=2:3;
doexp(1)=~all(size(tranii.img)==[501,501,164]);
try
doexp(2)=~all(size(cornii.img)==[501,501,164]);
catch
    doexp(2)=0;
end
doexp=exp(doexp);


% export images that have wrong dimension with correct bounding box
for export=doexp
    switch export
        case 2
            cf=        [options.root,options.prefs.patientdir,filesep,options.prefs.tranii,',1'];

            outf=options.prefs.tranii;
        case 3
            cf=        [options.root,options.prefs.patientdir,filesep,options.prefs.cornii,',1'];
            outf=options.prefs.cornii;
    end
    matlabbatch{1}.spm.util.imcalc.input = {[options.earoot,'templates',filesep,'bb.nii,1'];
        cf};
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
end
