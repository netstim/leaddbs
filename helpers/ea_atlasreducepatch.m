function atlases=ea_atlasreducepatch(atlases,red)

for atl=1:length(atlases.names)
    for side=1:2
        atlases.roi{atl,side}.fv=reducepatch(atlases.roi{atl,side}.fv,red);

    end
end

