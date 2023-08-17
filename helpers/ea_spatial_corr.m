function [R,p] = ea_spatial_corr(niftis,method,smooth)
% Correlates spatial correlations among a series of nifti files after
% conforming space to the first file.
if ~exist('method','var')
    method='Pearson';
end
if ~exist('smooth','var')
    smooth=0;
end


tmp=ea_getleadtempdir;
tmpniftis=niftis;
scnt=1;
R=nan(length(niftis),length(niftis),length(smooth));
p=R;
for s=smooth
    for n=1:length(niftis)
        uid=ea_generate_uuid;
        tmpniftis{n}=fullfile(tmp,[uid,'.nii']);
        copyfile(niftis{n},tmpniftis{n});
        if n>1
            ea_conformspaceto(tmpniftis{1},tmpniftis{n},1);
        end
        if s
            spm_smooth(tmpniftis{n},fullfile(tmp,['s',uid,'.nii']),[s s s]);
            nii=ea_load_nii(fullfile(tmp,['s',uid,'.nii']));
        else
            nii=ea_load_nii(tmpniftis{n});
            X(:,n)=nii.img(:);
        end
    end
    [R(:,:,scnt),p(:,:,scnt)]=corr(X,'Type',method,'rows','pairwise');
    scnt=scnt+1;
end

R=ea_nanmean(R,3);
p=ea_nanmean(p,3);