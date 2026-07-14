function ea_unified_draw_tracts(obj,vals,group,side,alphaind,fibcell,fibcmap,cmapind,usedidx)
% Plot fibers if any survived
if ~iscell(obj.drawvals)
    obj.drawvals = {};
end

if ~isempty(fibcell{group,side})
    prefs = ea_prefs;
    obj.drawobject.fiberfiltering{group,side} = streamtube(fibcell{group,side}, prefs.d3.fiberwidth);
    nones=repmat({'none'},size(fibcell{group,side}));
    [obj.drawobject.fiberfiltering{group,side}.EdgeColor]=nones{:};

    % Calulate fiber colors alpha values
    fibcolor = mat2cell(fibcmap{group}(cmapind{side},:), ones(size(fibcell{group,side})));
    fibalpha = mat2cell(alphaind{side},ones(size(fibcell{group,side})));

    % Set fiber colors and alphas
    [obj.drawobject.fiberfiltering{group,side}.FaceColor]=fibcolor{:};
    [obj.drawobject.fiberfiltering{group,side}.FaceAlpha]=fibalpha{:};
    obj.drawvals{group,side} = vals{group,side};
    obj.fiberdrawn.fibcell = fibcell;
    obj.fiberdrawn.vals = vals;
    obj.fiberdrawn.usedidx = usedidx;
else
    obj.drawobject.fiberfiltering{group,side} = {};
    obj.drawvals{group,side} = {};
end
end