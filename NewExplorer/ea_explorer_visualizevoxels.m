function [fibcell,vals,usedidx]=ea_explorer_visualizevoxels(obj)
%% Thresholding Part
if obj.thresholding.showsignificantonly
    obj.recentmodel.voxels.sigvals=ea_explorer_corrsignan(obj.recentmodel.voxels.vals,obj.recentmodel.voxels.pvals,obj);
    vals = obj.recentmodel.voxels.sigvals;
else
    vals = obj.recentmodel.voxels.vals;
end
for side=1:numel(vals)
    % Remove vals and fibers outside the thresholding range
    obj.stats.voxels.pos.available(side)=sum(cat(1,vals{:,side})>0);
    obj.stats.voxels.neg.available(side)=sum(cat(1,vals{:,side})<0);

    usedidx{1,side}=find(~isnan(vals{1,side})); % find nonnans
    vals{1,side}=vals{1,side}(usedidx{1,side}); % remove nans

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
    usedidx{1,side}(remove)=[];
    obj.stats.voxels.pos.shown(side)=sum(vals{1,side}>0);
    obj.stats.voxels.neg.shown(side)=sum(vals{1,side}<0);
end
%% Visualization Part
set(0,'CurrentFigure',obj.resultfig);

% reset colorbar
obj.colorbar.voxels=[];

allvals = full(vertcat(vals{1,:}));

if isempty(allvals) || all(isnan(allvals))
    % ea_cprintf('CmdWinWarnings', 'Empty or all-nan value found!\n');
    return;
end

if obj.thresholding.posvisible && all(allvals<0)
    obj.thresholding.posvisible = 0;
    fprintf('\n')
    warning('off', 'backtrace');
    warning('No positive voxels found, posvisible is set to 0 now!');
    warning('on', 'backtrace');
    fprintf('\n')
end

if obj.thresholding.negvisible && all(allvals>0)
    obj.thresholding.negvisible = 0;
    fprintf('\n')
    warning('off', 'backtrace');
    warning('No negative voxels found, negvisible is set to 0 now!');
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

if obj.thresholding.posvisible
    cmap = ea_colorgradient(gradientLevel, [1,1,1], obj.thresholding.poscolor);
    voxcmap.pos = ea_colorgradient(gradientLevel, cmap(shiftedCmapStart,:), obj.thresholding.poscolor);
    cmapind.pos = round(normalize(allvals,'range',[1,gradientLevel]));
    cmapind.pos=mat2cell(cmapind.pos, cellfun(@numel,vals))';
    % alphaind = normalize(1./(1+exp(-allvals)), 'range');
end
if obj.thresholding.negvisible
    cmap = ea_colorgradient(gradientLevel, obj.thresholding.negcolor, [1,1,1]);
    voxcmap.neg = ea_colorgradient(gradientLevel, obj.thresholding.negcolor, cmap(shiftedCmapEnd,:));
    cmapind.neg = round(normalize(allvals,'range',[1,gradientLevel]));
    cmapind.neg = mat2cell(cmapind.neg, cellfun(@numel,vals))';
    % alphaind = normalize(-1./(1+exp(-allvals)), 'range');
end

setappdata(obj.resultfig, ['voxcmap',obj.ID], voxcmap);


for side=1:size(vals,2)

    % Plot voxels if any survived
    %% Positive balls
    if any(vals{1,side} > 0) && obj.thresholding.posvisible
        posidx = usedidx{1,side}(vals{1,side}>0);
        posspot.nii=obj.results.space{1,side};
        posspot.nii.img(posidx)=vals{1,side}(vals{1,side}>0);
        posspot.name=['Positive' num2str(side)];
        posspot.niftiFilename=[posspot.name '.nii'];;
        posspot.binary=0;
        posspot.usesolidcolor=0;
        posspot.color=obj.thresholding.poscolor;
        posspot.colormap=voxcmap.pos;
        posspot.smooth=0;
        posspot.hullsimplify=0;
        posspot.threshold=0;
        obj.drawnsweetspots.pos{1,side}=ea_roi(posspot.niftiFilename,posspot);
    end
    if any(vals{1,side} < 0) && obj.thresholding.negvisible
        negidx = usedidx{1,side}(vals{1,side}<0);
        negspot.nii=obj.results.space{1,side};
        negspot.nii.img(negidx)=vals{1,side}(vals{1,side}<0);
        negspot.name=['Negative' num2str(side)];
        negspot.niftiFilename=[negspot.name '.nii'];
        negspot.binary=0;
        negspot.usesolidcolor=0;
        negspot.color=obj.thresholding.negcolor;
        negspot.colormap=voxcmap.neg;
        negspot.smooth=0;
        negspot.hullsimplify=0;
        negspot.threshold=0;
        obj.drawnsweetspots.neg{1,side}=ea_roi(negspot.niftiFilename,negspot);
    end
end


% Set colorbar tick positions and labels
% if ~isempty(allvals)
%     if obj.thresholding.posvisible && obj.thresholding.negvisible
%         tick = [1, floor(length(voxcmap)/2-40), ceil(length(voxcmap)/2+40), length(voxcmap)];
%         poscbvals = sort(allvals(allvals>0));
%         negcbvals = sort(allvals(allvals<0));
%         ticklabel = [negcbvals(1), negcbvals(end), poscbvals(1), poscbvals(end)];
%         ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
%     elseif obj.thresholding.posvisible
%         tick = [1, length(voxcmap)];
%         posvals = sort(allvals(allvals>0));
%         ticklabel = [posvals(1), posvals(end)];
%         ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
%     elseif obj.thresholding.negvisible
%         tick = [1, length(voxcmap)];
%         negvals = sort(allvals(allvals<0));
%         ticklabel = [negvals(1), negvals(end)];
%         ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
%     end
% end
% % store colorbar in object
% if exist('voxcmap','var') % could be no fibers present at all.
%     obj.colorbar.voxels.cmap = voxcmap;
%     obj.colorbar.voxels.tick = tick;
%     obj.colorbar.voxels.ticklabel = ticklabel;
% end
end