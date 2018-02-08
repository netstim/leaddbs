function sels=ea_storeupdatecortex(jtree,h)

for root=1:length(h.sg)
    for branch=1:length(h.sgsub) %side
        sels.branches{branch}=h.sgsub{branch}.getSelectionState;
        for leaf=1:length(h.sgsubfi{branch})
            sels.leaves{branch}{leaf}=h.sgsubfi{branch}{leaf}.getSelectionState;
            if isfield(h,'sgsubside')
                for side=1:length(h.sgsubside{branch}{leaf})
                    sels.sides{branch}{leaf}{side}=h.sgsubside{branch}{leaf}{side}.getSelectionState;
                end
            end
        end
    end
end
setappdata(jtree,'selectionstate',sels);