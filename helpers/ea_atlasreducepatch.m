function atlases=ea_atlasreducepatch(atlases,red)

for atl=1:length(atlases.names)
    for side=1:2
        atlases.fv{atl,side}=reducepatch(atlases.fv{atl,side},red);
        
    end
end

