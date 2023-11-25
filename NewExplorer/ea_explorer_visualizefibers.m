function [fibcell,vals,usedidx]=ea_explorer_visualizefibers(obj)
%% Thresholding Part
fibcell=obj.results.(ea_conn2connid(obj.connectome)).fibcell;

if obj.thresholding.showsignificantonly
    obj.recentmodel.fibers.sigvals=ea_explorer_corrsignan(obj.recentmodel.fibers.vals,obj.recentmodel.fibers.pvals,obj);
    vals = obj.recentmodel.fibers.sigvals;
else
    vals = obj.recentmodel.fibers.vals;
end
for side=1:numel(vals)
    % Remove vals and fibers outside the thresholding range
    obj.stats.fibers.pos.available(side)=sum(cat(1,vals{:,side})>0);
    obj.stats.fibers.neg.available(side)=sum(cat(1,vals{:,side})<0);

    usedidx{1,side}=find(~isnan(vals{1,side})); % find nonnans
    vals{1,side}=vals{1,side}(usedidx{1,side}); % remove nans
    fibcell{1,side}=fibcell{1,side}(usedidx{1,side}); % remove nans

    % separating and sort positive and negative vals
    switch obj.thresholding.threshstrategy
        case 'Fixed Amount' % here we want to create threshs for each side separately.
            posvals = sort(vals{1,side}(vals{1,side}>0),'descend');
            negvals = sort(vals{1,side}(vals{1,side}<0),'ascend');
        otherwise % in other cases, we want to apply the same thresh to both sides.
            allvals = vertcat(vals{1,:});
            posvals = sort(allvals(allvals>0),'descend');
            negvals = sort(allvals(allvals<0),'ascend');
            clear allvals
    end
    % positive thresholds
    if ~obj.thresholding.posvisible || ~obj.thresholding.showposamount(side) || isempty(posvals)
        posthresh = inf;
    else
        posthresh = ea_explorer_applythresholding(obj.thresholding.threshstrategy, posvals, obj.thresholding.showposamount(side));
    end
    % negative thresholds
    if ~obj.thresholding.negvisible || ~obj.thresholding.shownegamount(side) || isempty(negvals)
        negthresh = -inf;
    else
        negthresh = ea_explorer_applythresholding(obj.thresholding.threshstrategy, negvals, obj.thresholding.shownegamount(side));
    end
    remove = logical(logical(vals{1,side}<posthresh) .* logical(vals{1,side}>negthresh));
    vals{1,side}(remove)=[];
    fibcell{1,side}(remove)=[];
    usedidx{1,side}(remove)=[];
    obj.stats.fibers.pos.shown(side)=sum(vals{1,side}>0);
    obj.stats.fibers.neg.shown(side)=sum(vals{1,side}<0);
end
%% Visualization Part
set(0,'CurrentFigure',obj.resultfig);

% reset colorbar
obj.colorbar.fibers=[];

allvals = full(vertcat(vals{1,:}));

if isempty(allvals) || all(isnan(allvals))
    % ea_cprintf('CmdWinWarnings', 'Empty or all-nan value found!\n');
    return;
end

if obj.thresholding.posvisible && all(allvals<0)
    obj.thresholding.posvisible = 0;
    fprintf('\n')
    warning('off', 'backtrace');
    warning('No positive fibers found, posvisible is set to 0 now!');
    warning('on', 'backtrace');
    fprintf('\n')
end

if obj.thresholding.negvisible && all(allvals>0)
    obj.thresholding.negvisible = 0;
    fprintf('\n')
    warning('off', 'backtrace');
    warning('No negative fibers found, negvisible is set to 0 now!');
    warning('on', 'backtrace');
    fprintf('\n')
end

% colormap(gray);
gradientLevel = 1024;
cmapShiftRatio = 0.5;
shiftedCmapStart = round(gradientLevel*cmapShiftRatio)+1;
shiftedCmapEnd = gradientLevel-round(gradientLevel*cmapShiftRatio);
shiftedCmapLeftEnd = gradientLevel/2-round(gradientLevel/2*cmapShiftRatio);
shiftedCmapRightStart = round(gradientLevel/2*cmapShiftRatio)+1;

alphaind = ones(size(allvals));
if obj.thresholding.posvisible && obj.thresholding.negvisible
    cmap = ea_colorgradient(gradientLevel/2, obj.thresholding.negcolor, [1,1,1]);
    cmapLeft = ea_colorgradient(gradientLevel/2, obj.thresholding.negcolor, cmap(shiftedCmapLeftEnd,:));
    cmap = ea_colorgradient(gradientLevel/2, [1,1,1], obj.thresholding.poscolor);
    cmapRight = ea_colorgradient(gradientLevel/2, cmap(shiftedCmapRightStart,:), obj.thresholding.poscolor);
    fibcmap = [cmapLeft;cmapRight];
    cmapind = ones(size(allvals))*gradientLevel/2;
    cmapind(allvals<0) = round(normalize(allvals(allvals<0),'range',[1,gradientLevel/2]));
    cmapind(allvals>0) = round(normalize(allvals(allvals>0),'range',[gradientLevel/2+1,gradientLevel]));
    % alphaind(allvals<0) = normalize(-1./(1+exp(-allvals(allvals<0))), 'range');
    % alphaind(allvals>0) = normalize(1./(1+exp(-allvals(allvals>0))), 'range');
elseif obj.thresholding.posvisible
    cmap = ea_colorgradient(gradientLevel, [1,1,1], obj.thresholding.poscolor);
    fibcmap = ea_colorgradient(gradientLevel, cmap(shiftedCmapStart,:), obj.thresholding.poscolor);
    cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
    % alphaind = normalize(1./(1+exp(-allvals)), 'range');
elseif obj.thresholding.negvisible
    cmap = ea_colorgradient(gradientLevel, obj.thresholding.negcolor, [1,1,1]);
    fibcmap = ea_colorgradient(gradientLevel, obj.thresholding.negcolor, cmap(shiftedCmapEnd,:));
    cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
    % alphaind = normalize(-1./(1+exp(-allvals)), 'range');
end

setappdata(obj.resultfig, ['fibcmap',obj.ID], fibcmap);

cmapind=mat2cell(cmapind, cellfun(@numel,vals))';
alphaind=mat2cell(alphaind, cellfun(@numel,vals))';

for side=1:size(vals,2)
    % Plot fibers if any survived
    if ~isempty(fibcell{1,side})
        prefs = ea_prefs;
        if obj.fastrender
            tic
            pointthreshold = 10;
            numpoints = cellfun('size',fibcell{1,side},1);
            tmpidx=find(numpoints > pointthreshold);
            steps=ones(size(numpoints));
            steps(tmpidx)=numpoints(tmpidx)./pointthreshold;
            for i=1:numel(fibcell{1,side})
                fibcell{1,side}{i}=fibcell{1,side}{i}(round(1:steps(i):numpoints(i)),:);
            end
            obj.drawnstreamlines{1,side} = streamtube(fibcell{1,side}, prefs.d3.fiberwidth);
            disp(['Fast Render took: ' num2str(round(toc,2)) 's.'])
        else
            tic
            obj.drawnstreamlines{1,side} = streamtube(fibcell{1,side}, prefs.d3.fiberwidth);
            disp(['Slow Render took: ' num2str(round(toc,2)) 's.'])
        end
        set(obj.drawnstreamlines{1,side},'EdgeColor','none')
        % Calulate fiber colors alpha values
        fibcolor = mat2cell(fibcmap(cmapind{side},:), ones(size(fibcell{1,side})));
        fibalpha = mat2cell(alphaind{side},ones(size(fibcell{1,side})));
        % Set fiber colors and alphas
        [obj.drawnstreamlines{1,side}.FaceColor]=fibcolor{:};
        [obj.drawnstreamlines{1,side}.FaceAlpha]=fibalpha{:};
        obj.drawvals{1,side} = vals{1,side};
    else
        obj.drawnstreamlines{1,side} = {};
        obj.drawvals{1,side} = {};
    end
end

% Set colorbar tick positions and labels
if ~isempty(allvals)
    if obj.thresholding.posvisible && obj.thresholding.negvisible
        tick = [1, floor(length(fibcmap)/2-40), ceil(length(fibcmap)/2+40), length(fibcmap)];
        poscbvals = sort(allvals(allvals>0));
        negcbvals = sort(allvals(allvals<0));
        ticklabel = [negcbvals(1), negcbvals(end), poscbvals(1), poscbvals(end)];
        ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
    elseif obj.thresholding.posvisible
        tick = [1, length(fibcmap)];
        posvals = sort(allvals(allvals>0));
        ticklabel = [posvals(1), posvals(end)];
        ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
    elseif obj.thresholding.negvisible
        tick = [1, length(fibcmap)];
        negvals = sort(allvals(allvals<0));
        ticklabel = [negvals(1), negvals(end)];
        ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
    end
end
% store colorbar in object
if exist('fibcmap','var') % could be no fibers present at all.
    obj.colorbar.fibers.cmap = fibcmap;
    obj.colorbar.fibers.tick = tick;
    obj.colorbar.fibers.ticklabel = ticklabel;
end
end