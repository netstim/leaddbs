function ea_surfice_heatmap_menu(~,~,handles,sides)


[fis,path]=uigetfile({'.nii'},'Select files to visualize...','Multiselect','on');

for fi=1:length(fis)
    
   fullfi=fileparts(path,fis{fi}); 
    keyboard
    
end