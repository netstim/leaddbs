function ea_surfice_heatmap_menu(~,~,handles,sides,manualcontrast)
% exports heatmaps from selected niftis via surfice.


[fis,path]=uigetfile({'.nii'},'Select files to visualize...','Multiselect','on');

if ~iscell(fis)
    fis={fis};
end

for fi=1:length(fis)
    
    fullfi{fi}=fullfile(path,fis{fi});
end

if manualcontrast
    threshs=ea_sfc_getautothresh(fullfi);
end

[do,threshs]=ea_sfc_setthreshs(threshs,fis);
if strcmp(do,'proceed')
     ea_surficeoverlay_lr(fullfi,[],sides);
end