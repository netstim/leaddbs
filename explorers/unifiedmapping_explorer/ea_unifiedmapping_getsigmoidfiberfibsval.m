function fibsval = ea_unifiedmapping_getsigmoidfiberfibsval(obj, connid)
% Return stored sigmoid peak fiber values or request recalculation.

[available, resultField] = ...
    ea_unifiedmapping_hassigmoidfiberfibsval(obj, connid);
if ~available
    message = [ ...
        'Stored sigmoid fiber activation values are not available. ', ...
        'Recalculate Fiber Filtering to create them.'];
    ea_warndlg(message);
    error('UnifiedMapping:MissingSigmoidFiberValues', '%s', message);
end

fibsval = ...
    obj.results.fiberfiltering.(connid).(resultField).fibsval;
end
