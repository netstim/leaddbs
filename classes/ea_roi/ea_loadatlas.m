function atlases=ea_loadatlas(atlname,resultfig,ht)

if ~exist('resultfig','var')
    resultfig=gcf;
end
if ~exist('ht','var')
    ht=[];
end

load([ea_space([],'atlases'),atlname,filesep,'atlas_index.mat']);
for atl=1:size(atlases.roi,1)
    for side=1:2
        if ~isa(atlases.roi{atl,side},'ea_roi')
            atlases.roi{atl,side}=ea_struct2roi(atlases.roi{atl,side},resultfig,ht);
        end
    end
end


