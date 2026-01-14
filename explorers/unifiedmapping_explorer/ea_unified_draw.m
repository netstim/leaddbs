function ea_unified_draw(obj,vals,fibcell,usedidx) %for cv live visualize
%function draw(obj,vals,fibcell)
% re-define plainconn (since we do not store it)
%this is all based on which tool is selected
%match draw tool with field name
if strcmp(obj.drawTool,'sweetspotmapping')
    flag2tool = 'sweetspotmapping';
    if isempty(obj.drawobject)
        obj.drawobject.sweetspotmapping = struct;
    end
elseif strcmp(obj.drawTool,'fiberfiltering')
    flag2tool = 'fiberfiltering';
    if isempty(obj.drawobject)
        obj.drawobject.fiberfiltering = struct;
    end
elseif strcmp(obj.drawTool,'networkmapping')
    flag2tool = 'networkmapping';
    if isempty(obj.drawobject)
        obj.drawobject.networkmapping = struct;
    end
end
[vals,fibcell,usedidx] = ea_unifiedmapping_calcstats(obj);
if strcmp(obj.drawTool,'fiberfiltering')
    ea_updatePAM_VAT(obj);

end

export=nan;
%ea_discfibers_showroi(obj);
obj.fiberdrawn.fibcell = fibcell;
obj.fiberdrawn.vals = vals;
obj.fiberdrawn.usedidx = usedidx;
% if show connected (white) fibers, also calculate those

allvals{1}=[]; % need to use a loop here - cat doesnt work in all cases with partly empty cells..
if size(vals,2)==2 % can be a single cell in case of custom code (pseudoM setting).
    allvals{2}=[];
end
for v=1:size(vals,1)
    allvals{1}=[allvals{1};vals{v,1}];
    if size(vals,2)==2
        allvals{2}=[allvals{2};vals{v,2}];
    end
end
obj.stats.pos.shown(1)=sum(allvals{1}>0);
obj.stats.neg.shown(1)=sum(allvals{1}<0);
if size(vals,2)>1 % bihemispheric usual case
    obj.stats.pos.shown(2)=sum(allvals{2}>0);
    obj.stats.neg.shown(2)=sum(allvals{2}<0);
end

set(0,'CurrentFigure',obj.resultfig);

domultitract=size(vals,1)>1; % if color by groups is set will be positive.
if ~isfield(obj.M,'groups')
    obj.M.groups.group=1;
    obj.M.groups.color=ea_color_wes('all');
end

switch obj.multitractmode
    case 'Split & Color By Group'
        linecols=obj.M.groups.color;
    case 'Split & Color By Subscore'
        linecols = obj.subscore.colors;
    case 'Split & Color By PCA'
        linecols = obj.subscore.pcacolors;
    case 'Single Tract Analysis'
        linecols = obj.M.groups.color;
end
if isfield(obj.drawobject,flag2tool)
    if isempty(obj.drawobject.(flag2tool)) % check if prior object has been stored
        obj.drawobject.(flag2tool)=getappdata(obj.resultfig,['dt_',obj.ID]); % store handle of tract to figure.
    end
    obj.drawobject.(flag2tool)={};
else
    obj.drawobject.(flag2tool) = {};
end


if strcmp(obj.multitractmode,'Single Tract Analysis') || obj.subscore.special_case
    % reset colorbar
    obj.colorbar=[];
    if ~any([obj.posvisible,obj.negvisible])
        return
    end
end

for group=1:size(vals,1) % vals will have 1x2 in case of bipolar drawing and Nx2 in case of group-based drawings (where only positives are shown).
    % vals will also be >1 for subscore tracts
    % Vertcat all values for colorbar construction
    if domultitract && ~obj.subscore.special_case
        if ~any([obj.subscore.posvisible(group),obj.subscore.negvisible(group)])
            continue
        end
    end
    if domultitract
        obj.subscore.vis.pos_shown(group,1)=sum(vals{group,1}>0);
        obj.subscore.vis.neg_shown(group,1)=sum(vals{group,1}<0);
        if size(vals{group},2)==1 % unihemispheric case, should then be pseudoM case
            obj.subscore.vis.pos_shown(group,2) = 0;
            obj.subscore.vis.neg_shown(group,2) = 0;
        elseif (size(vals{group,2},1))>1 % bihemispheric usual case
            obj.subscore.vis.pos_shown(group,2)=sum(vals{group,2}>0);
            obj.subscore.vis.neg_shown(group,2)=sum(vals{group,2}<0);
        elseif length(vals(group,:)) == 2 % in the case that it still exist
            obj.subscore.vis.pos_shown(group,2) = 0;
            obj.subscore.vis.neg_shown(group,2) = 0;
        end
    end
    allvals = full(vertcat(vals{group,:}));

    if isempty(allvals) || all(isnan(allvals))
        % ea_cprintf('CmdWinWarnings', 'Empty or all-nan value found!\n');
        continue;
    else
        allvals(isnan(allvals)) = 0;
        allvals(isinf(allvals)) = 0; % ignore infs for colormap generation.
    end

    if strcmp(obj.multitractmode,'Split & Color By Subscore') || strcmp(obj.multitractmode,'Split & Color By PCA')
        if obj.subscore.special_case && group == 1
            vals_across_subscores = full(vertcat(vals{1:size(vals,1),:}));
            check_and_update_visibility(obj, 'posvisible', vals_across_subscores, '<0');
            check_and_update_visibility(obj, 'negvisible', vals_across_subscores, '>0');
        elseif obj.subscore.special_case
            check_and_update_visibility(obj.subscore, 'posvisible', allvals, '<0', group);
            check_and_update_visibility(obj.subscore, 'negvisible', allvals, '>0', group);
        else
            check_and_update_visibility(obj, 'posvisible', allvals, '<0');
            check_and_update_visibility(obj, 'negvisible', allvals, '>0');
        end
    end
    colormap(gray);
    gradientLevel = 1024;
    cmapShiftRatio = 0.5;
    shiftedCmapStart = round(gradientLevel*cmapShiftRatio)+1;
    shiftedCmapEnd = gradientLevel-round(gradientLevel*cmapShiftRatio);
    shiftedCmapLeftEnd = gradientLevel/2-round(gradientLevel/2*cmapShiftRatio);
    shiftedCmapRightStart = round(gradientLevel/2*cmapShiftRatio)+1;

    if obj.subscore.special_case || ~domultitract
        pos_active = obj.posvisible;
        neg_active = obj.negvisible;
    else
        pos_active = obj.subscore.posvisible(group);
        neg_active = obj.subscore.negvisible(group);
    end
    if domultitract
        switch obj.multitractmode
            %logic is different for groups (pos & neg cannot be
            %shown together), whereas for PCA it is not as
            %such. Therefore, I am splitting the cases into
            %two.
            case 'Split & Color By Group'
                obj.poscolor = linecols(group,:);
                obj.negcolor = [0.94,0.97,1.00];
            case 'Split & Color By Subscore'
                obj.poscolor = obj.subscore.colors{1,2}(group,:); % positive main color
                obj.negcolor = obj.subscore.colors{1,1}(group,:); % negative main color
            case 'Split & Color By PCA'
                obj.poscolor = obj.subscore.pcacolors(group,:);
                obj.negcolor = [0.94,0.97,1.00];
        end
    else
        obj.poscolor = obj.poscolor; % positive main color
        obj.negcolor = obj.negcolor; % negative main color
    end
    cmap = ea_colorgradient(gradientLevel, [1,1,1], obj.poscolor);
    if pos_active && ~neg_active
        fibcmap{group} = ea_colorgradient(gradientLevel, cmap(shiftedCmapStart,:), obj.poscolor);
        cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
        alphaind = ones(size(allvals));
        % alphaind = normalize(allvals, 'range');
    elseif ~pos_active && neg_active
        cmap = ea_colorgradient(gradientLevel, obj.negcolor, [1,1,1]);
        fibcmap{group} = ea_colorgradient(gradientLevel, obj.negcolor, cmap(shiftedCmapEnd,:));
        cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
        alphaind = ones(size(allvals));
        % alphaind = normalize(-allvals, 'range');
    elseif pos_active && neg_active
        if strcmp(obj.multitractmode,'Split & Color By Group')
            warndlg(sprintf(['Please choose either "Show Positive Fibers" or "Show Negative Fibers".',...
                '\nShow both positive and negative fibers is not supported when "Color by Subscore Variable" is on.']));
            return;
        else
            cmap = ea_colorgradient(gradientLevel/2, obj.negcolor, [1,1,1]);
            cmapLeft = ea_colorgradient(gradientLevel/2, obj.negcolor, cmap(shiftedCmapLeftEnd,:));
            cmap = ea_colorgradient(gradientLevel/2, [1,1,1], obj.poscolor);
            cmapRight = ea_colorgradient(gradientLevel/2, cmap(shiftedCmapRightStart,:), obj.poscolor);
            fibcmap{group} = [cmapLeft;cmapRight];
            cmapind = ones(size(allvals))*gradientLevel/2;
            cmapind(allvals<0) = round(normalize(allvals(allvals<0),'range',[1,gradientLevel/2]));
            cmapind(allvals>0) = round(normalize(allvals(allvals>0),'range',[gradientLevel/2+1,gradientLevel]));
            alphaind = ones(size(allvals));
        end
    end
    setappdata(obj.resultfig, ['fibcmap',obj.ID], fibcmap);
    if size(vals,2)>1 % standard case
        cmapind = mat2cell(cmapind, [numel(vals{group,1}), numel(vals{group,2})])';
        alphaind = mat2cell(alphaind, [numel(vals{group,1}), numel(vals{group,2})])';
    else % potential scripting case, only one side
        cmapind = mat2cell(cmapind, numel(vals{group,1}))';
        alphaind = mat2cell(alphaind, numel(vals{group,1}))';
    end
    for side=1:size(vals,2)
        switch obj.drawTool
            case {'sweetspotmapping','networkmapping'}
                ea_unified_draw_volumetric(obj,vals,group,side,gradientLevel,fibcmap);
            case 'fiberfiltering'
                ea_unified_draw_tracts(obj,vals,group,side,alphaind,fibcell,fibcmap,cmapind);
                %  obj.activate_tractset; % function to highlight tracts ac
        end
        if domultitract % introduce small jitter for visualization
            fibcell{group,side}=ea_unifiedmapping_addjitter(fibcell{group,side},0.03);
        end
    end
    % Set colorbar tick positions and labels
    if ~isempty(allvals)
        if pos_active && neg_active
            tick{group} = [1, floor(length(fibcmap{group})/2-40), ceil(length(fibcmap{group})/2+40), length(fibcmap{group})];
            poscbvals = sort(allvals(allvals>0));
            negcbvals = sort(allvals(allvals<0));
            if ~isempty(negcbvals) && ~isempty(poscbvals)
                ticklabel{group} = [negcbvals(1), negcbvals(end), poscbvals(1), poscbvals(end)];
                ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
            else
                continue
            end
        elseif pos_active && ~neg_active
            tick{group} = [1, length(fibcmap{group})];
            posvals = sort(allvals(allvals>0));
            ticklabel{group} = [posvals(1), posvals(end)];
            ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
        elseif neg_active && ~pos_active
            tick{group} = [1, length(fibcmap{group})];
            negvals = sort(allvals(allvals<0));
            ticklabel{group} = [negvals(1), negvals(end)];
            ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
        end
        % store colorbar in object
        if exist('fibcmap','var') % could be no fibers present at all.
            obj.colorbar.cmap = fibcmap;
            obj.colorbar.tick = tick;
            obj.colorbar.ticklabel = ticklabel;
        end

    end

end
    function copyobj(thisObj,newObj)
        % Construct a new object based on a deep copy of the current
        % object of this class by copying properties over.
        props = properties(thisObj);
        for i = 1:length(props)
            % Use Dynamic Expressions to copy the required property.
            % For more info on usage of Dynamic Expressions, refer to
            % the section "Creating Field Names Dynamically" in:
            % web([docroot '/techdoc/matlab_prog/br04bw6-38.html#br1v5a9-1'])
            newObj.(props{i}) = thisObj.(props{i});
        end
    end

end






function ea_updatePAM_VAT(obj)
% Extract connectome ID once
conn_id = ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome);

% Ensure 'plainconn' has correct connectivity values
try
    key = "PAM_Ttest";
    if obj.calcsettings.connectivity_type ~= 2
        key = "VAT_Ttest";
    end
    obj.results.fiberfiltering.(conn_id).('plainconn').fibsval = obj.results.fiberfiltering.(conn_id).(key).fibsval;
catch
    warn_msg = "Connectivity indices were not stored. Please recalculate or stay with the same model (VAT or PAM)";
    ea_warndlg(warn_msg);
    disp(repmat("=", 1, 55)); disp("WARNING: " + warn_msg); disp(repmat("=", 1, 55));
end

% If switching connectivity, update fiber pathways
if obj.calcsettings.switch_connectivity == 1
    cfile = find_connectivity_file(obj);
    load(cfile, 'fibers', 'idx');

    try
        for side = 1:2
            conn_key = "connFiberInd_PAM";
            if obj.calcsettings.connectivity_type ~= 2
                conn_key = "connFiberInd_VAT";
            end

            connFiberInd = obj.results.fiberfiltering.(conn_id).(conn_key){side};
            connFiber = fibers(ismember(fibers(:,4), connFiberInd), 1:3);
            obj.results.fiberfiltering.(conn_id).fibcell{side} = mat2cell(connFiber, idx(connFiberInd));
        end
    catch
        ea_warndlg(warn_msg);
        disp(repmat("=", 1, 55)); disp("WARNING: " + warn_msg); disp(repmat("=", 1, 55));
    end
end

% Ensure 'totalFibers' is stored
if ~isfield(obj.results.fiberfiltering.(conn_id), 'totalFibers')
    cfile = find_connectivity_file(obj);
    load(cfile, 'idx');
    obj.results.fiberfiltering.(conn_id).totalFibers = length(idx);
end

% Helper function to find the connectivity file
    function cfile = find_connectivity_file(obj)
        if obj.calcsettings.multi_pathways == 1
            cfile = fullfile(fileparts(obj.analysispath), obj.calcsettings.fibfilt_connectome, 'merged_pathways.mat');
        else
            cfile = fullfile(ea_getconnectomebase('dMRI'), obj.calcsettings.fibfilt_connectome, 'data.mat');
        end
    end
end