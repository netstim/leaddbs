feats={'ct','bg','bs','cb'};

for f=1:length(feats)
    %copyfile(['atropos_',feats{f},'.nii'],[feats{f},'.nii']);
    %ea_conformspaceto('222.nii',[feats{f},'.nii']);
    nii=ea_load_nii([feats{f},'.nii']);
    idx{f}=find(nii.img(:)>0.2);
    
end

save('feature_idx','idx','feats');


