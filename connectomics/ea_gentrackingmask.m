function ea_gentrackingmask(options,threshold)
directory=[options.root,options.patientname,filesep];

ea_newseg(directory,options.prefs.prenii_unnormalized,0,options);

%% Coreg options.prefs.prenii_unnormalized to b0 (for label.mat and FTR-Normalization)
copyfile([directory,options.prefs.prenii_unnormalized],[directory,'c',options.prefs.prenii_unnormalized]);
copyfile([directory,'c2',options.prefs.prenii_unnormalized],[directory,'cc2',options.prefs.prenii_unnormalized]);
ea_coreg2images(options,[directory,'c',options.prefs.prenii_unnormalized],[directory,options.prefs.b0],[directory,'c',options.prefs.prenii_unnormalized],{[directory,'cc2',options.prefs.prenii_unnormalized]},1);

movefile([directory,'cc2',options.prefs.prenii_unnormalized],[directory,'trackingmask.nii']);
delete([directory,'c',options.prefs.prenii_unnormalized]);

tr=ea_load_nii([options.root,options.patientname,filesep,'trackingmask.nii']);
if threshold
    tr.img=tr.img>0.8;
    tr.fname=[options.root,options.patientname,filesep,'ttrackingmask.nii'];
    ea_write_nii(tr);
end