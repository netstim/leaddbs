function atlases=ea_reorderatlas(atlases,order)
% reorders an atlas struct based on the order supplied
try atlases.names=atlases.names(order); end
try atlases.types=atlases.types(order); end
try atlases.tissuetypes=atlases.tissuetypes(order); end
atlases.colors=atlases.colors(order);
try atlases.fv=atlases.fv(order,:); end
try atlases.roi=atlases.roi(order,:); end

try atlases.cdat=atlases.cdat(order,:); end
try atlases.XYZ=atlases.XYZ(order,:); end
try atlases.pixdim=atlases.pixdim(order,:); end
try atlases.colorc=atlases.colorc(order,:); end
try atlases.normals=atlases.normals(order,:); end

try
for p=1:length(atlases.presets)
    atlases.presets(p).show=find(ismember(order,atlases.presets(p).show));
    atlases.presets(p).hide=find(ismember(order,atlases.presets(p).hide));
end
end
try
for l=1:length(atlases.labels)
atlases.labels{l}=atlases.labels{l}(order);
end
end