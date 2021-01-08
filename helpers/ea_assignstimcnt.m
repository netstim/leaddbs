function [ea_stats,thisstim]=ea_assignstimcnt(ea_stats,S)
thisstim=1;

if isfield(ea_stats,'stimulation')
    if ~isfield(ea_stats.stimulation,'label') % legacy, delete older stimulations
        ea_stats=rmfield(ea_stats,'stimulation');
    else
        for stimfield=1:length(ea_stats.stimulation)
           if strcmp(ea_stats.stimulation(stimfield).label,S.label)
              thisstim=stimfield;
              return
           end
           % if not assigned yet, assign here.
           thisstim=length(ea_stats.stimulation)+1;
        end
    end
end

ea_stats.stimulation(thisstim).label=S.label;
