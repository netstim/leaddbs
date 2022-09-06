function peak = ea_discfibers_getpeak(vals, posvisible, negvisible)
% Get peak fiber values
% ea_discfibers_getpeak(vals, obj.posvisible, obj.negvisible)

if ~exist('posvisible', 'var')
    posvisible = 1;
end

if ~exist('negvisible', 'var')
    negvisible = 1;
end

if posvisible && ~negvisible
    peak = ea_nanmax(vals, 1);
elseif ~posvisible && negvisible
    peak = ea_nanmin(vals, 1);
elseif posvisible && negvisible
    pospeak = ea_nanmax(vals, 1);
    negpeak = ea_nanmin(vals, 1);
    peak = ea_nansum([pospeak; negpeak]);
end
