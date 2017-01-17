function ea_genwarp2space(fromwhich)

if ~exist('fromwhich','var')
    fromwhich='MNI_ICBM_2009b_NLIN_ASYM';
end
if strcmp(ea_getspace,fromwhich)
    return
end
load([ea_space,'ea_space_def.mat'])

mkdir([ea_space,'subj']);
for t=1:length(spacedef.templates)
    copyfile([ea_space,spacedef.templates{t},'.nii'], [ea_space,'subj',filesep,'anat_',spacedef.templates{t},'.nii']);
end


originalspace=ea_getspace;


% do normalization of that 'subject':
options.root=ea_space;
options.patientname='subj';
options.prefs=ea_prefs('');
options.earoot=ea_getearoot;
% temporarily switch back to from space to get warps
ea_switchspace([],[],fromwhich,'mute');
[options,presentfiles]=ea_assignpretra(options);
ea_normalize_spmshoot(options);


% delete all anat files:
delete([ea_space,'subj',filesep,'anat_*.nii']);
movefile([ea_space,'subj',filesep,'y_ea_normparams.nii'],[ea_space,'y_ea_normparams.nii']);
movefile([ea_space,'subj',filesep,'y_ea_inv_normparams.nii'],[ea_space,'y_ea_inv_normparams.nii']);
rmdir([ea_space,'subj'],'s');

% switch space back to the user selected one

ea_switchspace([],[],originalspace,'mute');
