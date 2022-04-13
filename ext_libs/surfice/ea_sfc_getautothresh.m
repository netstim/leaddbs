function threshold = ea_sfc_getautothresh(fis)
% Calculate thresholds for visualization in SufIce

threshold = nan(length(fis), 4);

for fi=1:length(fis)
    nii = ea_load_nii(fis{fi});
    nii.img(nii.img==0) = nan;
    imgstd = ea_nanstd(nii.img(:));
    if min(nii.img(:))<0 % has negative values, calculate 4 vals
        threshold(fi,1) = 0.5*imgstd;
        threshold(fi,2) = 2.7*imgstd;
        threshold(fi,3) = -2.7*imgstd;
        threshold(fi,4) = -0.5*imgstd;
    else % only positive vals
        threshold(fi,1) = 0.04*imgstd;
        threshold(fi,2) = 1.9*imgstd;
        threshold(fi,3) = nan;
        threshold(fi,4) = nan;
    end

    if ~any(nii.img(:)>0)
        threshold(fi, 1:2) = nan;
    end
end
