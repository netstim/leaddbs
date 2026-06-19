function [AllX] = ea_unified_networkmapping_recalcvals_sk(obj, sk, basefield)

if ~exist('basefield', 'var') || isempty(basefield)
    basefield = 'connval';
end

AllX=obj.results.networkmapping.(ea_conn2connid(obj.calcsettings.netmap_connectome)).(basefield);

td=ea_getleadtempdir;
% get nifti space
switch size(AllX,2)
    case 902629 % 222
        space=ea_load_nii(fullfile(ea_getearoot,'templates','spacedefinitions','222.nii.gz'));
    case 7221032 % 111
        space=ea_load_nii(fullfile(ea_getearoot,'templates','spacedefinitions','111.nii.gz'));
    case 69402312 % 555
        space=ea_load_nii(fullfile(ea_getearoot,'templates','spacedefinitions','555.nii.gz'));
end
%% Load in nifti files as matrix
for s=1:size(AllX,1)
    space.img(:)=AllX(s,:);
    uid=ea_generate_uuid;
    space.fname=fullfile(td,[uid,'.nii']);
    ea_write_nii(space);
    [X]=ea_genX({space.fname},[],[],ea_unified_getobjmask(obj),sk);    
    AllX(s,:) = X;
end
