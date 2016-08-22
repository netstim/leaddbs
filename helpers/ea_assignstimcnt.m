function [ea_stats,thisstim]=ea_assignstimcnt(ea_stats,S)
if isfield(ea_stats,'stimulation')
    if ~isfield(ea_stats.stimulation,'label') % legacy, delete older stimulations
        ea_stats=rmfield(ea_stats,'stimulation');
        thisstim=1;
    else
        for stimfield=1:length(ea_stats.stimulation)
           if strcmp(ea_stats.stimulation(stimfield).label,S.label)
              thisstim=stimfield;
              return
           end
        end
    end
else % no stim has been done at all before
    thisstim=1;
end

ea_stats.stimulation(thisstim).label=S.label;