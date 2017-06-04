function [ixs,ixt]=ea_getsubindex(sel,sidec,surfs,togglebuttons)

sel = char(sel);

if regexp(sel, '<HTML><BODY>')
    [sel ~] = regexp(sel,'(?:&nbsp;)+(\S*)</FONT></BODY></HTML>','tokens','match');
    sel = sel{1}{1};
end

for p=1:length(surfs)
    if strcmp(surfs(p).Tag,[sel,sidec])
        ixs=p;
        break
    end
end

for p=1:length(surfs)
    if strcmp(togglebuttons(p).Tag,[sel,sidec])
        ixt=p;
        break
    end
end