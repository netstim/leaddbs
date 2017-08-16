function atlases=ea_removeatlascomponents(atlases,toremove)
% removes an atlas entry based on the ix supplied


reorder=1:length(atlases.names);
reorder(toremove)=[];
atlases=ea_reorderatlas(atlases,reorder);
