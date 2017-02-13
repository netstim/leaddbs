function ea_genwarp2space(fromwhich)

if strcmp(ea_getspace,fromwhich)
    return
end
if exist([ea_getspace,filesep,fromwhich],'dir') % warp exists
    return
end

foreignspaceroot=[ea_getearoot,'templates',filesep,'space',filesep,fromwhich,filesep];
load([foreignspaceroot,'ea_space_def.mat'])

mkdir([ea_space,fromwhich]);
for t=1:length(spacedef.templates)
    copyfile([foreignspaceroot,spacedef.templates{t},'.nii'], [ea_space,fromwhich,filesep,'anat_',spacedef.templates{t},'.nii']);
end




% do normalization of that 'subject':
options.root=ea_space;
options.patientname=fromwhich;
options.prefs=ea_prefs('');
options.earoot=ea_getearoot;
options.modality=1;
% temporarily switch back to from space to get warps
[options,presentfiles]=ea_assignpretra(options);

% perform three runs to refine ? this needs to be really precise and is worth the additional time!
ea_normalize_spmdartel(options);
ea_normalize_spmdartel(options);
ea_normalize_spmdartel(options);


