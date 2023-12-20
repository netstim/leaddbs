
%% Cleartune Parameters
obj.cleartune.results = abj.cleartuneresults;% copy of results for auto tuning functions
obj.cleartune.efields = abj.cleartuneefields;% efields used to calc results
obj.cleartune.injected = abj.cleartuneinjected;% status to report file has injected values
obj.cleartune.optim = abj.CleartuneOptim;
obj.cleartune.vars = abj.cleartunevars;
%% Visualization
obj.visualization.roivisible = abj.roivisible; % show ROI (usually VTAs)
obj.visualization.connfibvisible = abj.connfibvisible;% show all connected tracts in white
obj.visualization.posvisible = abj.posvisible; % pos tract visible
obj.visualization.negvisible = abj.negvisible; % neg tract visible
obj.visualization.showposamount = abj.showposamount; % two entries for right and left
obj.visualization.shownegamount = abj.shownegamount; % two entries for right and left
obj.visualization.poscolor = abj.poscolor; % positive main color
obj.visualization.negcolor = abj.negcolor; % negative main color
obj.visualization.showsignificantonly = abj.showsignificantonly;
%% Multipath
obj.multipath.multi_pathways = abj.multi_pathways ; % if structural connectome is devided into pathways (multiple .mat in dMRI_MultiTract)
obj.multipath.map_list = abj.map_list;% list that contains global indices of the first fibers in each pathway (relevant when multi_pathways = 1)
obj.multipath.pathway_list = abj.pathway_list; % list that contains names of pathways (relevant when multi_pathways = 1