function atlases=ea_loadatlas(atlname,resultfig,ht)
% Load atlases and do struct2roi conversion when needed

load([ea_space([],'atlases'),atlname,filesep,'atlas_index.mat'], 'atlases');

if isfield(atlases, 'roi')
    if ~exist('resultfig','var')
        resultfig=gcf;
    end

    if ~exist('ht','var')
        ht=[];
    end

    for atl=1:size(atlases.roi,1)
        for side=1:size(atlases.roi,2)
            if ~isa(atlases.roi{atl,side},'ea_roi')
                atlases.roi{atl,side}=ea_struct2roi(atlases.roi{atl,side},resultfig,ht);
            end
        end
    end
end
