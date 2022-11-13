function [vals,fibcell,usedidx] = ea_discfibers_loadModel_calcstats(obj,vals_connected)

% vals_connected are defined in the subspace of connected fiber for the
% particular lead-group loaded in Fiber Filtering Explorer

% vals does not contain zero entries
vals = cell(size(vals_connected));
usedidx = cell(size(vals_connected));
for voter = 1:size(vals,1)
    for side=1:size(vals,2)

        usedidx{voter,side}=find(vals_connected{voter,side});
        vals{voter,side}=vals_connected{voter,side}(usedidx{voter,side}); % final weights for surviving fibers

        %fibcell{group,side}=obj.results.(ea_conn2connid(obj.connectome)).fibcell{side}(~isnan(vals{group,side}));
        fibcell{voter,side}=obj.results.(ea_conn2connid(obj.connectome)).fibcell{side}(usedidx{voter,side});
    end
end


