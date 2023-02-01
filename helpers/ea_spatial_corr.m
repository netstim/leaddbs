function [R,p] = ea_spatial_corr(niftis,method)
% Correlates spatial correlations among a series of nifti files after
% conforming space to the first file.
if ~exist('method','var')
    method='Pearson';
end

tmp=ea_getleadtempdir;
tmpniftis=niftis;
for n=1:length(niftis)
    uid=ea_generate_uuid;
    tmpniftis{n}=fullfile(tmp,[uid,'.nii']);
    copyfile(niftis{n},tmpniftis{n});
    if n>1
        ea_conformspaceto(tmpniftis{1},tmpniftis{n},1);
    end
    nii=ea_load_nii(tmpniftis{n});
    X(:,n)=nii.img(:);
end

[R,p]=corr(X,'Type',method,'rows','pairwise');
