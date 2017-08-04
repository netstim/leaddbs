function atlases=ea_reorderatlas(atlases,order)
% reorders an atlas struct based on the order supplied
atlases.names=atlases.names(order);
atlases.types=atlases.types(order);
atlases.tissuetypes=atlases.tissuetypes(order);
atlases.colors=atlases.colors(order);
atlases.fv=atlases.fv(order,:);
atlases.cdat=atlases.cdat(order,:);
atlases.XYZ=atlases.XYZ(order,:);
atlases.pixdim=atlases.pixdim(order,:);
atlases.colorc=atlases.colorc(order,:);
atlases.normals=atlases.normals(order,:);


for p=1:length(atlases.presets)
    atlases.presets(p).show=find(ismember(order,atlases.presets(p).show));
    atlases.presets(p).hide=find(ismember(order,atlases.presets(p).hide));
end
for l=1:length(atlases.labels)
atlases.labels{l}=atlases.labels{l}(order);
end