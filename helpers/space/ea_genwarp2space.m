function ea_genwarp2space(fromwhich)

if strcmp(ea_getspace,fromwhich)
    return
end
foreignspaceroot=[ea_getearoot,'templates',filesep,'space',filesep,fromwhich,filesep];
load([foreignspaceroot,'ea_space_def.mat'])

mkdir([ea_space,'subj']);
for t=1:length(spacedef.templates)
    copyfile([foreignspaceroot,spacedef.templates{t},'.nii'], [ea_space,'subj',filesep,'anat_',spacedef.templates{t},'.nii']);
end




% do normalization of that 'subject':
options.root=ea_space;
options.patientname='subj';
options.prefs=ea_prefs('');
options.earoot=ea_getearoot;
options.modality=1;
% temporarily switch back to from space to get warps
[options,presentfiles]=ea_assignpretra(options);

ea_normalize_spmshoot(options);


