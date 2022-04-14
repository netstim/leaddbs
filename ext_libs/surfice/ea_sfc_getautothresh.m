function threshold = ea_sfc_getautothresh(images)
% Calculate thresholds for visualization in SufIce

if ischar(images)
    images = {images};
end

threshold = nan(length(images), 4);

for f=1:length(images)
    nii = ea_load_nii(images{f});
    nii.img(nii.img==0) = nan;
    imgstd = ea_nanstd(nii.img(:));
    if min(nii.img(:))<0 % has negative values, calculate 4 vals
        threshold(f,1) = 0.5*imgstd;
        threshold(f,2) = 2.7*imgstd;
        threshold(f,3) = -0.5*imgstd;
        threshold(f,4) = -2.7*imgstd;
    else % only positive vals
        threshold(f,1) = 0.04*imgstd;
        threshold(f,2) = 1.9*imgstd;
        threshold(f,3) = nan;
        threshold(f,4) = nan;
    end

    if ~any(nii.img(:)>0)
        threshold(f, 1:2) = nan;
    end
end
