function ea_genwarp2space(fromwhich)

if strcmp(ea_getspace,fromwhich)
    return
end

if exist([ea_space,fromwhich],'dir') && exist([ea_space,fromwhich,filesep,'glanatComposite',ea_getantstransformext([ea_space,fromwhich,filesep])],'file') % warp exists
    return
end

foreignspaceroot=[ea_getearoot,'templates',filesep,'space',filesep,fromwhich,filesep];
load([foreignspaceroot,'spacedef.mat'])

mkdir([ea_space,fromwhich]);
for t=1:length(spacedef.templates)
    if exist([foreignspaceroot,spacedef.templates{t},'.nii'],'file')
    copyfile([foreignspaceroot,spacedef.templates{t},'.nii'], [ea_space,fromwhich,filesep,'anat_',spacedef.templates{t},'.nii']);
    end
end




% do normalization of that 'subject':
options.root=ea_space;
options.patientname=fromwhich;
options.prefs=ea_prefs('');
options.earoot=ea_getearoot;
options.modality=1;
options.coregmr.method='ANTs (Avants 2008)';
% temporarily switch back to from space to get warps
[options,presentfiles]=ea_assignpretra(options);

norm_method_applied{1}='ea_normalize_ants';
save([ea_space,fromwhich,filesep,'ea_normmethod_applied.mat'],'norm_method_applied');


% % perform three runs to refine ? this needs to be really precise and is worth the additional time!
% ea_normalize_spmdartel(options);
% ea_normalize_spmdartel(options);
% ea_normalize_spmdartel(options);
ea_normalize_ants(options);


