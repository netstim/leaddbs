function [vals,fibcell,usedidx] = ea_unifiedmapping_calcstats(obj,patsel,Iperm)

if ismethod(obj, 'repair_loaded_explorer')
    obj.repair_loaded_explorer;
end

%No fibcell and usedidx for sweetspotmapping or networkmapping
switch obj.drawTool
    case {'sweetspotmapping','networkmapping'}
        fibcell = {};
        usedidx = {};
end
% NB: for PCA, we are going to reassign I later in the function
if ~exist('Iperm','var')
    I=obj.responsevar;
else % used in permutation based statistics - in this case the real improvement can be substituted with permuted variables.
    I=Iperm;
end

% quickly recalc stats
if ~exist('patsel','var') % patsel can be supplied directly (in this case, obj.patientselection is ignored), e.g. for cross-validations.
    patsel=obj.patientselection;
end
if size(I,2)==1 % 1 entry per patient, not per electrode
    I=[I,I]; % both sides the same;
end

if obj.mirrorsides
    I=[I;I];
end

if ~isempty(obj.covars)
    for i=1:length(obj.covars)
        if obj.mirrorsides
            covars{i} = [obj.covars{i};obj.covars{i}];
        else
            covars{i} = obj.covars{i};
        end
        for side=1:2
            if size(covars{i},2)==1 % single entry, use the same for both sides
                I(:,side)=ea_resid(covars{i},I(:,side));
            elseif size(covars{i},2)==2 % two entries in covars, use each for the respective side:
                I(:,side)=ea_resid(covars{i}(:,side),I(:,side));
            end
        end
    end  
end
groups=1;
dogroups = 0;
dosubscores = 0;
switch obj.drawTool
    case 'sweetspotmapping' %sweetspotmapping
        init_val = obj.results.sweetspotmapping.efield';
        for i=1:length(init_val)
            init_val{i} = init_val{i}';
        end

    case 'fiberfiltering' %fiberfiltering
        connid = ea_conn2connid(obj.calcsettings.fibfilt_connectome);        
        if obj.calcsettings.connectivity_type == 2
            selectedFibcell = ...
                obj.results.fiberfiltering.(connid).pam_fibers.fibcell;
        
        elseif strcmp(obj.e_field_metric, 'Projection')
            selectedFibcell = ...
                obj.results.fiberfiltering.(connid).efield_proj.fibcell;
        
        else
            selectedFibcell = ...
                obj.results.fiberfiltering.(connid).efield_fibers.fibcell;
        end

        if obj.calcsettings.connectivity_type == 2
            switch obj.statsettings.stimulationmodel
                case 'VTA'
                    init_val = obj.results.fiberfiltering.(connid).('PAM_Ttest').fibsval;
                otherwise
                    init_val = obj.results.fiberfiltering.(connid).('PAM_probA').fibsval;
            end
                
        else
            switch obj.statsettings.stimulationmodel
                case 'Sigmoid Field'
                    fibsval_raw = obj.results.fiberfiltering.(connid).(ea_unifiedmapping_method2methodid(obj)).fibsval;
                    init_val = fibsval_raw;  % initialize
                    for side = 1:size(fibsval_raw,2)
                        init_val{1,side}(:,:) = ea_unified_probabilityActivationFunction(fibsval_raw{1,side}(:,:));
                    end
                otherwise
                    init_val = cellfun(@full, obj.results.fiberfiltering.(connid).(ea_unifiedmapping_method2methodid(obj)).fibsval, 'Uni', 0);
            end
        end

    case 'networkmapping' %networkmapping
        init_val = ea_get_AllX(obj);
        init_val = {init_val'}; %side will always be 1 - both the hemispheres are combined for this analysis!
end

for group = groups
    gval = init_val;
    gpatsel=patsel;
    if obj.mirrorsides
        gpatsel=[gpatsel,gpatsel+length(obj.allpatients)];
    end
    if size(I,2)==1 % 1 entry per patient, not per electrode
        I=[I,I]; % both sides the same;
    end
    if obj.mirrorsides
        I=[I;I];
    end
    for side=1:numel(gval)
        vals{group,side}=nan(size(gval{side}(:,gpatsel),1),1);
        if obj.showsignificantonly
            pvals{group,side}=vals{group,side};
        end
        switch obj.drawTool
            case 'sweetspotmapping'
                switch obj.statsettings.stimulationmodel
                    case 'VTA'
                        % Sweet-spot calculation stores the thresholded input
                        % maps for all stimulation models. Convert retained
                        % voxels to a binary VTA here so two-sample tests see
                        % explicit connected (1) and unconnected (0) groups.
                        vtaValues = gval{side}(:,gpatsel);
                        missingValues = isnan(vtaValues);
                        vtaValues = double(vtaValues > 0);
                        vtaValues(missingValues) = nan;
                        gval{side}(:,gpatsel) = vtaValues;
                        Nmap=ea_nansum(gval{side}(:,gpatsel),2);
                        gval{side}(Nmap<((obj.statsettings.connthreshold/100)*length(gpatsel)),gpatsel)=nan;
                    case 'Electric Field'
                        % Optional legacy behavior for reproducing analyses
                        % that treated low E-field values as missing. The
                        % current default keeps these values (including zero)
                        % unless the user explicitly sets a NaN threshold.
                        if isfield(obj.statsettings, 'nanthreshold') && ...
                                ~isempty(obj.statsettings.nanthreshold)
                            gval{side}(gval{side} <= obj.statsettings.nanthreshold) = nan;
                        end
                        Nmap=ea_nansum((gval{side}(:,gpatsel)>obj.statsettings.efieldthreshold_spot),2);
                        gval{side}(Nmap<round((obj.statsettings.connthreshold/100)*length(gpatsel)),gpatsel)=nan;
                    case 'Sigmoid Field'
                        gval{side}(:, gpatsel) = ea_SigmoidFromEfield(gval{side}(:, gpatsel));
                        Nmap = ea_nansum( gval{side}(:, gpatsel) > obj.statsettings.efieldthreshold_spot, 2);
                        minN = round((obj.statsettings.connthreshold / 100) * length(gpatsel)); 
                        gval{side}(Nmap < minN, gpatsel) = nan;

                end
                %initialize vals and pvals if necessary

                %rules for removing nonempty values, since there are many nans in
                %the voxel wise analysis we can choose two different ways of
                %dealing with the vals

                % Two-sample tests require both connected and unconnected
                % observations. Apply the upper coverage threshold before
                % constructing valsin so overly common voxels are excluded from
                % the data passed to the statistical test.
                if strcmpi(obj.statsettings.statfamily, '2-Sample Tests')
                    maxConnected = (1 - obj.statsettings.connthreshold/100) * length(gpatsel);
                    gval{side}(Nmap > maxConnected, gpatsel) = nan;
                end

                nonempty = sum(gval{side}(:,gpatsel),2,'omitnan')>0;
                nonemptyidx=find(nonempty);
                valsin=gval{side}(nonempty,gpatsel);
            case 'fiberfiltering'
                    switch obj.statsettings.stimulationmodel
                    case 'VTA'
                        Nmap=ea_nansum(gval{side}(:,gpatsel),2);
                        gval{side}(Nmap<((obj.statsettings.connthreshold/100)*length(gpatsel)),gpatsel)=nan;
                    case 'Sigmoid Field'
                        pafThreshold = obj.statsettings.efieldthreshold_tract;
                        Nmap=ea_nansum((gval{side}(:,gpatsel)>pafThreshold),2);
                        gval{side}(Nmap < round((obj.statsettings.connthreshold/100) * length(gpatsel)), gpatsel) = nan;
                    case 'Electric Field'
                        Nmap=ea_nansum((gval{side}(:,gpatsel)>obj.statsettings.efieldthreshold_tract),2);
                        gval{side}(Nmap<round((obj.statsettings.connthreshold/100)*length(gpatsel)),gpatsel)=nan;
                    end
                
                %initialize vals and pvals if necessary

                %rules for removing nonempty values, since there are many nans in
                %the voxel wise analysis we can choose two different ways of
                %dealing with the vals

                % Two-sample tests require both connected and unconnected
                % observations. Apply the upper coverage threshold before
                % constructing valsin so overly common fibers are excluded from
                % the data passed to the statistical test.
                if strcmpi(obj.statsettings.statfamily, '2-Sample Tests')
                    maxConnected = (1 - obj.statsettings.connthreshold/100) * length(gpatsel);
                    gval{side}(Nmap > maxConnected, gpatsel) = nan;
                end

                nonempty = sum(gval{side}(:,gpatsel),2,'omitnan')>0;
                nonemptyidx=find(nonempty);
                valsin=gval{side}(nonempty,gpatsel);
            case 'networkmapping'
              
                allVals = gval{side}(:, gpatsel);  % extract the submatrix
                % change the % from GUI to r
                corrThreshold = obj.statsettings.efieldthreshold_network / 100;

                % minimum patients required
                minN = ceil((obj.statsettings.connthreshold / 100) * size(allVals,2));
                        
                
                NmapPos = sum(allVals >= corrThreshold, 2, 'omitnan');
                NmapNeg = sum(allVals <= -corrThreshold, 2, 'omitnan');
               
                keep = NmapPos >= minN | NmapNeg >= minN;
                
                % remove voxels not meeting the criteria
                allVals(~keep, :) = nan;
                gval{side}(:, gpatsel) = allVals;
                valsin = allVals;                 
        end

        
        %define outcome
        outcomein=I(gpatsel,side);

        %start the stats process
        if strcmp(obj.statsettings.statfamily,'Correlations')
            disp([obj.drawTool,': Calculating ' obj.statsettings.stattest ' correlation for side ' num2str(side) '...']);
        else
            
            disp(['Calculating ' obj.statsettings.stattest ' for side ' num2str(side) '...']);
        end
        %define all the stattests available
        stattests=ea_explorer_statlist(obj.responsevar);
        %determine which is the test you are going to run
        [is,idx] = ismember(obj.statsettings.stattest, stattests.name);
        if ~is
            ea_error(['Function for test ',obj.statsettings.statset,' missing.']);
        end

        % ranksum removes NaN outcomes and fails when either the connected
        % or unconnected group is empty. Exclude untestable features here so
        % the shared stat-test implementation does not need to be changed.
        if strcmp(stattests.file(idx), 'ea_explorer_stats_ranksumtest') && ...
                ismember(obj.drawTool, {'sweetspotmapping','fiberfiltering'})
            finiteOutcome = isfinite(outcomein(:)).';
            hasConnectedGroup = any((valsin == 1) & finiteOutcome, 2);
            hasUnconnectedGroup = any((valsin == 0) & finiteOutcome, 2);
            testable = hasConnectedGroup & hasUnconnectedGroup;

            valsin = valsin(testable, :);
            nonemptyidx = nonemptyidx(testable);
        end

        %do the actual calculation
        if ~isempty(valsin) && ~isempty(outcomein)
            [valsout,psout]=feval(stattests.file(idx),valsin,outcomein,obj.statsettings.H0); % apply test
            switch  obj.drawTool
                case {'sweetspotmapping','fiberfiltering'}
                    vals{group,side}(nonemptyidx)=valsout;
                    if exist('pvals','var')
                        pvals{group,side}(nonemptyidx)=psout;
                    end
                case 'networkmapping'
                    vals{group,side}=valsout;
                    if exist('pvals','var')
                        pvals{group,side}=psout;
                    end
            end
          
        else
            if isempty(valsin)
                ea_cprintf('CmdWinWarnings', 'group %d side %d: empty valsin!\n', group, side);
            end
            if isempty(outcomein)
                ea_cprintf('CmdWinWarnings', 'group %d side %d: empty outcomein!\n', group, side);
            end
        end
    end
end

if obj.showsignificantonly
    vals=ea_corrsignan(vals,pvals,obj);
end

if strcmp(obj.threshstrategy,'Unthresholded')
    % Unthresholded means no magnitude/rank threshold, but the selected
    % sign visibility still defines which values belong to the model. This
    % mask is also used by cross-validation through the returned vals.
    for group = groups
        for side = 1:numel(gval)
            if dosubscores || dogroups
                if obj.subscore.special_case
                    showPositive = obj.posvisible;
                    showNegative = obj.negvisible;
                else
                    showPositive = obj.subscore.posvisible(group);
                    showNegative = obj.subscore.negvisible(group);
                end
            else
                showPositive = obj.posvisible;
                showNegative = obj.negvisible;
            end

            if ~showPositive
                vals{group,side}(vals{group,side} > 0) = nan;
            end
            if ~showNegative
                vals{group,side}(vals{group,side} < 0) = nan;
            end
        end
    end

    if strcmp(obj.drawTool,'fiberfiltering')
        for side=1:numel(gval)
            usedidx{group,side} = find(isfinite(vals{group,side}));
            fibcell{group,side} = selectedFibcell{side}(usedidx{group,side});            
            vals{group,side} = vals{group,side}(usedidx{group,side}); % final weights for surviving fibers
            if exist('pvals','var')
                pvals{group,side} = pvals{group,side}(usedidx{group,side}); % final weights for surviving fibers
            end

            obj.stats.pos.available(side)=sum(cat(1,vals{:,side})>0);
            obj.stats.neg.available(side)=sum(cat(1,vals{:,side})<0);

            if dosubscores || dogroups
                if ~obj.subscore.special_case
                    obj.subscore.vis.pos_available(group,side)=sum(cat(1,vals{group,side})>0);
                    obj.subscore.vis.neg_available(group,side)=sum(cat(1,vals{group,side})<0);
                end
            end
        end
    end
   return
else
    switch obj.drawTool
        case {'sweetspotmapping','networkmapping'}
        
        for group = groups
            for side=1:numel(gval)
                obj.stats.pos.available(side)=sum(cat(1,vals{:,side})>0);
                obj.stats.neg.available(side)=sum(cat(1,vals{:,side})<0);
                if dosubscores || dogroups
                    if ~obj.subscore.special_case
                        obj.subscore.vis.pos_available(group,side)=sum(cat(1,vals{group,side})>0);
                        obj.subscore.vis.neg_available(group,side)=sum(cat(1,vals{group,side})<0);
                    end
                end
            end
        end
        case 'fiberfiltering'
        % Clean up non-finite values from fibcell and vals
        for group=groups
            for side=1:numel(gval)
                usedidx{group,side} = find(isfinite(vals{group,side}));
                fibcell{group,side} = selectedFibcell{side}(usedidx{group,side});                vals{group,side} = vals{group,side}(usedidx{group,side}); % final weights for surviving fibers
                if exist('pvals','var')
                    pvals{group,side} = pvals{group,side}(usedidx{group,side}); % final weights for surviving fibers
                end

                obj.stats.pos.available(side)=sum(cat(1,vals{:,side})>0);
                obj.stats.neg.available(side)=sum(cat(1,vals{:,side})<0);

                if dosubscores || dogroups
                    if ~obj.subscore.special_case
                        obj.subscore.vis.pos_available(group,side)=sum(cat(1,vals{group,side})>0);
                        obj.subscore.vis.neg_available(group,side)=sum(cat(1,vals{group,side})<0);
                    end
                end
            end
        end
    end
    unthresholdedVals = vals; % Need to keep this original vals when calculating the (same) threshold to both sides.

    % Thresholding
    for group=groups
        for side=1:numel(gval)
            switch obj.threshstrategy
                case 'Fixed Amount' % here we want to create thresholds for each side separately.
                    posvals = sort(unthresholdedVals{group,side}(unthresholdedVals{group,side}>0),'descend');
                    negvals = sort(unthresholdedVals{group,side}(unthresholdedVals{group,side}<0),'ascend');

                otherwise % in other cases, we want to apply the same threshold to both sides.
                    allvals = vertcat(unthresholdedVals{group,:});
                    posvals = sort(allvals(allvals>0),'descend');
                    negvals = sort(allvals(allvals<0),'ascend');
            end

            % positive thresholds
            if dosubscores || dogroups
                if obj.subscore.special_case
                    if ~obj.posvisible || ~obj.showposamount(side) || isempty(posvals)
                        posthresh = inf;
                    else
                        posthresh = ea_fibValThresh(obj.threshstrategy, posvals, obj.showposamount(side));
                    end
                else
                    if ~obj.subscore.posvisible(group) || ~obj.subscore.vis.showposamount(group,side) || isempty(posvals)
                        posthresh = inf;
                    else
                        posthresh = ea_fibValThresh(obj.threshstrategy, posvals, obj.subscore.vis.showposamount(group,side));
                    end
                end
            else
                if ~obj.posvisible || ~obj.showposamount(side) || isempty(posvals)
                    posthresh = inf;
                else
                    posthresh = ea_fibValThresh(obj.threshstrategy, posvals, obj.showposamount(side));
                end
            end

            % negative thresholds
            if dosubscores || dogroups
                if obj.subscore.special_case
                    if ~obj.negvisible || ~obj.shownegamount(side) || isempty(negvals)
                        negthresh = -inf;
                    else
                        negthresh = ea_fibValThresh(obj.threshstrategy, negvals, obj.shownegamount(side));
                    end
                else
                    if ~obj.subscore.negvisible(group) || ~obj.subscore.vis.shownegamount(group,side) || isempty(negvals)
                        negthresh = -inf;
                    else
                        negthresh = ea_fibValThresh(obj.threshstrategy, negvals, obj.subscore.vis.shownegamount(group,side));
                    end
                end
            else
                if ~obj.negvisible || ~obj.shownegamount(side) || isempty(negvals)
                    negthresh = -inf;
                else
                    negthresh = ea_fibValThresh(obj.threshstrategy, negvals, obj.shownegamount(side));
                end
            end

            if ~obj.runwhite 
                % Remove vals and fibers outside the thresholding range (set by
                % sliders)
                remove = vals{group,side}<posthresh & vals{group,side}>negthresh;
                vals{group,side}(remove)=nan;
                if strcmp(obj.drawTool,'fiberfiltering')       
                    vals{group,side}(remove)=[];
                    fibcell{group,side}(remove)=[];
                    usedidx{group,side}(remove)=[];
                end
            end
        end
    end

end
function AllX=ea_get_AllX(obj)
connid = ea_conn2connid(obj.calcsettings.netmap_connectome);
basefield = 'connval';
cacheRoot = obj.results.networkmapping.(connid);
if isprop(obj, 'mask_vta_fp') && obj.mask_vta_fp
    if isfield(cacheRoot, 'connval_vtamasked')
        basefield = 'connval_vtamasked';
        cachefield = 'vtamasked';
    else
        warning('VTA-masked network fingerprints were requested, but connval_vtamasked is missing. Falling back to unmasked fingerprints.');
        cachefield = '';
    end
else
    cachefield = '';
end

addchar='';
if obj.smooth_fp
    addchar=[addchar,'s'];
end
if obj.normalize_fp
    addchar=[addchar,'k'];
end
if isempty(addchar) % no s, no k
    AllX=cacheRoot.(basefield);
    return
end
try
    if isempty(cachefield)
        AllX=obj.results.networkmapping.(connid).(ea_conn2connid(lower(obj.cvmask))).(addchar).connval;
    else
        AllX=obj.results.networkmapping.(connid).(cachefield).(ea_conn2connid(lower(obj.cvmask))).(addchar).connval;
    end
catch
    [AllX] = ea_unified_networkmapping_recalcvals_sk(obj, addchar, basefield);
    if isempty(cachefield)
        obj.results.networkmapping.(connid).(ea_conn2connid(lower(obj.cvmask))).(addchar).connval=AllX;
    else
        obj.results.networkmapping.(connid).(cachefield).(ea_conn2connid(lower(obj.cvmask))).(addchar).connval=AllX;
    end
end

function fibValThreshold = ea_fibValThresh(threshstrategy, vals, threshold)
    switch threshstrategy
        case 'Percentage Relative to Peak'
            range = vals(1) - vals(end);
            fibValThreshold = vals(1) - threshold/100 * range;
            if range == 0
                if vals(1) > 0
                    fibValThreshold = fibValThreshold - eps*10;
                else
                    fibValThreshold = fibValThreshold + eps*10;
                end
            end
        case 'Percentage Relative to Amount'
            index = round((threshold/100)*length(vals));
            if index <=0
                fibValThreshold = vals(1);
            else
                fibValThreshold = vals(index);
            end
        case 'Fixed Amount'
            if length(vals)>round(threshold)
                fibValThreshold=vals(round(threshold));
            else
                fibValThreshold=vals(end);
            end
        case 'Histogram (CDF)'
            if vals(1) > 0
                [fx, x] = ecdf(vals);
                fibValThreshold = x(find(fx>=(1-threshold), 1));
            else
                [fx, x] = ecdf(-vals);
                fibValThreshold = -x(find(fx>=(1-threshold), 1));
            end
        case 'Fixed Fiber Value'
            fibValThreshold = threshold;
    end



function result = ea_isnan(input_array,flag)
    if size(input_array,2) > 1
        if strcmp(flag,'any')
            op = any(isnan(input_array));
        elseif strcmp(flag,'all')
            op = all(isnan(input_array));
        end
        if length(find(op)) > 1
            result = 1;
        else
            result = 0;

        end
    else
        if strcmp(flag,'any')
            result = any(isnan(input_array));
        elseif strcmp(flag,'all')
            result = all(isnan(input_array));
        end
    end

return
