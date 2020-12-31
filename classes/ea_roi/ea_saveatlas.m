function ea_saveatlas(atlname,atlases)

if isfield(atlases, 'roi')
    disp(atlases.roi)
     for roi=1:size(atlases.roi,1)
        for side=1:size(atlases.roi,2)
            atlases.roi{roi,side}=ea_roi2struct(atlases.roi{roi,side});
        end
    end
end

save([ea_space([],'atlases'),atlname,filesep,'atlas_index.mat'],'atlases','-v7.3');
