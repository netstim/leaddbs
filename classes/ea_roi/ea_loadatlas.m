function atlases = ea_loadatlas(atlases,resultfig,ht)
% Load atlases and do struct2roi conversion when needed

if ischar(atlases) % Input is not atlases struct itself
    if isfile(atlases) % Input is atlas_index.mat file path
        load(atlases);
    else % Input is atlas name (suppose in MNI space)
        load([ea_space([],'atlases'),atlases,filesep,'atlas_index.mat'], 'atlases');
    end
end

if isfield(atlases, 'roi')
    if isstruct(atlases.roi)
        atlases.roi = arrayfun(@(x) {x}, atlases.roi);
    end

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
