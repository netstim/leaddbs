function [isbinary,minmax]=ea_discfibers_checkcustomNii(vatlist)
% small function to check whether input custom ROI (pseudoM case) are
% binary and/or return their min/max values

for v=1:length(vatlist)

    nii=ea_load_nii(vatlist{v});

    isbinary(v)=ea_isbinary(nii.img);
    minmax(v,:)=[ea_nanmin(nii.img(:)), ea_nanmax(nii.img(:))];
end

isbinary=all(isbinary);
minmax=[ea_nanmin(minmax(:,1)),ea_nanmax(minmax(:,2))];





