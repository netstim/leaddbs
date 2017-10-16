function ea_tonemapct_file(options,nativemni)
ea_dispt('Tonemapping CT...');
if ~exist('nativemni','var')
    nativemni='native';
else
    if isempty(nativemni)
            nativemni='native';
    end
end
switch nativemni
    case 'native'
        directory=[options.root,options.patientname,filesep];
        if exist([directory,options.prefs.ctnii_coregistered],'file')
            ct=ea_load_nii([directory,options.prefs.ctnii_coregistered]);
            ct.fname=[directory,'tp_',options.prefs.ctnii_coregistered];
            ct.img=ea_tonemap_ct(ct.img);
            ea_write_nii(ct);
        else
            fprintf('%s not present. Skipping tonemapping.\n', options.prefs.ctnii_coregistered);
        end
    case 'mni'
        directory=[options.root,options.patientname,filesep];
        if exist([directory,options.prefs.gctnii],'file')
            ct=ea_load_nii([directory,options.prefs.gctnii]);
            ct.fname=[directory,'tp_',options.prefs.gctnii];
            ct.img=ea_tonemap_ct(ct.img);
            ea_write_nii(ct);
        else
            fprintf('%s not present. Skipping tonemapping.\n', options.prefs.gctnii);
        end
end

ea_dispt('');