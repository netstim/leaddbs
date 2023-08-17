function peak = ea_discfibers_getpeak(vals, posvisible, negvisible, peakmode)
% Get peak fiber values
% ea_discfibers_getpeak(vals, obj.posvisible, obj.negvisible, 'peak')

if ~exist('posvisible', 'var')
    posvisible = 1;
end

if ~exist('negvisible', 'var')
    negvisible = 0;
end

if ~exist('peakmode', 'var')
    peakmode = 'peak';
end

switch peakmode
    case {1, 'peak'}
        endInd = 1;
    case {0.05, 'peak5'}
        endInd = ceil(size(vals,1).*0.05);
    otherwise
        if isnumeric(peakmode) && peakmode<1
            endInd = ceil(size(vals,1).*peakmode);
        else
            error('Invalid peak mode!');
        end
end

if posvisible && ~negvisible
    vals = sort(vals, 'descend', 'MissingPlacement', 'last');
    vals(vals<0) = nan; % Remove negative values
    peak = ea_nansum(vals(1:endInd,:), 1);
elseif ~posvisible && negvisible
    vals = sort(vals, 'ascend');
    vals(vals>0) = nan; % Remove positive values
    peak = ea_nansum(vals(1:endInd,:), 1);
elseif posvisible && negvisible
    posvals = sort(vals, 'descend', 'MissingPlacement', 'last');
    posvals(posvals<0) = nan; % Remove negative values
    negvals = sort(vals, 'ascend');
    negvals(negvals>0) = nan; % Remove negative values
    peak = ea_nansum([posvals(1:endInd,:); negvals(1:endInd,:)], 1);
else
    error('Both positive and negative values are deselected!');
end
