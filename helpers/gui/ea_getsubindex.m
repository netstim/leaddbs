function [ixs,ixt]=ea_getsubindex(sel,sidec,surfs,togglebuttons,uselabelname,atlases)

sel = char(sel);

if regexp(sel, '<HTML><BODY>')
    [sel, ~] = regexp(sel,'(?:&nbsp;)+(.*)</FONT></BODY></HTML>','tokens','match');
    sel = sel{1}{1};
end
if uselabelname ~= 0
    [~, nsel] = ismember(sel,atlases.labels{uselabelname});
    if isfield(atlases,'tissuetypes') && (atlases.tissuetypes(nsel)==2)
        % Use fileparts to extract the name, fiber atlas files are
        % always named as XXX.mat
        [~, sel] = fileparts(atlases.names{nsel});
    else
        [~, sel] = ea_niifileparts(atlases.names{nsel});
    end
end

for p=1:length(surfs)
    if strcmp(surfs{p}.Tag,[sel,sidec])
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