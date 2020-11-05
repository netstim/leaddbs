function ea_saveatlas(atlname,atlases)


for roi=1:size(atlases.roi,1)
    for side=1:2
        atlases.roi{roi,side}=ea_roi2struct(atlases.roi{roi,side});
    end
end



save([ea_space([],'atlases'),atlname,filesep,'atlas_index.mat'],'atlases','-v7.3');
