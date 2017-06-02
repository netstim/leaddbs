function [ixs,ixt]=ea_getsubindex(sel,sidec,surfs,togglebuttons)

for p=1:length(surfs)
    if strcmp(surfs(p).Tag,[char(sel),sidec])
        ixs=p;
        break
    end
end

for p=1:length(surfs)
    if strcmp(togglebuttons(p).Tag,[char(sel),sidec])
        ixt=p;
        break
    end
end