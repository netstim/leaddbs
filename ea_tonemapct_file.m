function ea_tonemapct_file(options)
ea_dispt('Tonemapping CT...');
directory=[options.root,options.patientname,filesep];
ct=ea_load_nii([directory,options.prefs.ctnii_coregistered]);
ct.fname=[directory,'tp_',options.prefs.ctnii_coregistered];

ct.img=ea_tonemap_ct(ct.img);
ea_write_nii(ct);

ea_dispt('');