function R=ea_spatialsimilarity(fixedim,testims,mask)


fixed=ea_load_nii(fixedim);
if exist('mask','var')
    msk=ea_load_nii(mask);
else
    msk=fixed;
    msk.img(:)=1;
end
mix=logical(msk.img(:));
X=fixed.img(mix);

for t=1:length(testims)
    test=ea_load_nii(testims{t});
    X=[X,test.img(mix)];
end

R=corr(X,'rows','pairwise');
R=R(2:end,1);
