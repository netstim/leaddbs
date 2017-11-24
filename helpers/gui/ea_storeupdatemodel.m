function sels=ea_storeupdatemodel(jtree,h)


for branch=1:length(h.sg)
    sels.branches{branch}=h.sg{branch}.getSelectionState;
    for leaf=1:length(h.sgsub{branch})
        sels.leaves{branch}{leaf}=h.sgsub{branch}{leaf}.getSelectionState;
        if isfield(h,'sgsubside')
            for side=1:length(h.sgsubside{branch}{leaf})
                sels.sides{branch}{leaf}{side}=h.sgsubside{branch}{leaf}{side}.getSelectionState;
            end
        end
    end
end

setappdata(jtree,'selectionstate',sels);