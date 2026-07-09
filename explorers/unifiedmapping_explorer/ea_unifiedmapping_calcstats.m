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
            init_val = obj.results.fiberfiltering.(connid).('PAM_probA').fibsval;
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
                        Nmap=ea_nansum(gval{side}(:,gpatsel),2);
                        gval{side}(Nmap<((obj.statsettings.connthreshold/100)*length(gpatsel)),gpatsel)=nan;
                    otherwise
                        % old method; only if variable in workspace exists
                        if evalin('base','exist(''threshold_method'',''var'')')
                            threshold_method = evalin('base','threshold_method');
                            switch threshold_method
                                case 'old_method' % as was implemented in lead dbs v3.2.1
                                    gval{side}(gval{side}<=obj.statsettings.nanthreshold) = nan;
                                    Nmap=ea_nansum((gval{side}(:,gpatsel)>obj.statsettings.efieldthreshold),2);
                                    gval{side}(Nmap<round((obj.statsettings.connthreshold/100)*length(gpatsel)),gpatsel)=nan;
                            end
                        else %here we use the new percentile method
                            % compute global percentile across all selected patients
                            allVals = gval{side}(:, gpatsel);  % extract the submatrix
                            allVals(allVals == 0) = NaN;
                            thr = prctile(allVals(:), obj.statsettings.efieldthreshold);  % flatten and compute percentile
                            % mask fibers above threshold
                            Nmap = ea_nansum(gval{side}(:,gpatsel) >= thr, 2);
                            % apply connection threshold
                            gval{side}(Nmap < round((obj.statsettings.connthreshold/100) * length(gpatsel)), gpatsel) = nan;
                        end  
                end
                %initialize vals and pvals if necessary

                %rules for removing nonempty values, since there are many nans in
                %the voxel wise analysis we can choose two different ways of
                %dealing with the vals

                nonempty = sum(gval{side}(:,gpatsel),2,'omitnan')>0;
                nonemptyidx=find(nonempty);
                valsin=gval{side}(nonempty,gpatsel);
            case 'fiberfiltering'
                    switch obj.statsettings.stimulationmodel
                    case 'VTA'
                        Nmap=ea_nansum(gval{side}(:,gpatsel),2);
                        gval{side}(Nmap<((obj.statsettings.connthreshold/100)*length(gpatsel)),gpatsel)=nan;
                    case 'Sigmoid Field'
                        pafThreshold = obj.statsettings.efieldthreshold;
                        Nmap=ea_nansum((gval{side}(:,gpatsel)>pafThreshold),2);
                        gval{side}(Nmap < round((obj.statsettings.connthreshold/100) * length(gpatsel)), gpatsel) = nan;
                    otherwise
                        % old method; only if variable in workspace exists
                        if evalin('base','exist(''threshold_method'',''var'')')
                            threshold_method = evalin('base','threshold_method');
                            switch threshold_method
                                case 'old_method' % as was implemented in lead dbs v3.2.1
                                    gval{side}(gval{side}<=obj.statsettings.nanthreshold) = nan;
                                    Nmap=ea_nansum((gval{side}(:,gpatsel)>obj.statsettings.efieldthreshold),2);
                                    gval{side}(Nmap<round((obj.statsettings.connthreshold/100)*length(gpatsel)),gpatsel)=nan;
                            end
                        else %here we use the new percentile method
                            % compute global percentile across all selected patients
                            allVals = gval{side}(:, gpatsel);  % extract the submatrix
                            allVals(allVals == 0) = NaN;
                            thr = prctile(allVals(:), obj.statsettings.efieldthreshold);  % flatten and compute percentile
                            % mask fibers above threshold
                            Nmap = ea_nansum(gval{side}(:,gpatsel) >= thr, 2);
                            % apply connection threshold
                            gval{side}(Nmap < round((obj.statsettings.connthreshold/100) * length(gpatsel)), gpatsel) = nan;
                        end  
                    end
                
                %initialize vals and pvals if necessary

                %rules for removing nonempty values, since there are many nans in
                %the voxel wise analysis we can choose two different ways of
                %dealing with the vals

                nonempty = sum(gval{side}(:,gpatsel),2,'omitnan')>0;
                nonemptyidx=find(nonempty);
                valsin=gval{side}(nonempty,gpatsel);
            case 'networkmapping'
                if evalin('base','exist(''threshold_method'',''var'')')
                    threshold_method = evalin('base','threshold_method');
                        switch threshold_method
                            case 'old_method'
                                valsin = gval{side}(:,gpatsel);
                        end

                else %here we use the old method
                        % % compute global percentile across all patients
                        allVals = gval{side}(:, gpatsel);  % extract the submatrix
                        % % allVals(allVals == 0) = NaN;
                        % thr = prctile(allVals(:), obj.statsettings.efieldthreshold);  % flatten and compute percentile
                        % % mask voxels above threshold
                        % Nmap = ea_nansum(gval{side}(:,gpatsel) >= thr, 2);
                        % gval{side}(Nmap < round((obj.statsettings.connthreshold/100) * length(gpatsel)), gpatsel) = nan;
                        % valsin = gval{side}(:, gpatsel);
                        % 
                        % allVals = gval{side}(:, gpatsel);

                        % Separate signed network values. Direction matters here, so do not
                        % threshold positives and negatives together.
                        posVals = allVals(allVals > 0);
                        negVals = allVals(allVals < 0);
                
                        posMask = false(size(allVals));
                        negMask = false(size(allVals));
                
                        if ~isempty(posVals)
                            posThr = prctile(posVals, obj.statsettings.efieldthreshold);
                            posMask = allVals >= posThr;
                        end
                
                        if ~isempty(negVals)
                            negThr = prctile(negVals, 100 - obj.statsettings.efieldthreshold);
                            negMask = allVals <= negThr;
                        end
                
                        NmapPos = ea_nansum(posMask, 2);
                        NmapNeg = ea_nansum(negMask, 2);
                
                        minN = round((obj.statsettings.connthreshold/100) * length(gpatsel));
                        keep = NmapPos >= minN | NmapNeg >= minN;
                
                        gval{side}(~keep, gpatsel) = nan;
                        valsin = gval{side}(:, gpatsel);
                end  
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
        %do the actual calculation
        if ~isempty(valsin) && ~isempty(outcomein)
            if strcmp(obj.statsettings.statfamily,'2-sample Tests')
                % only in case of VTAs (given two-sample-t-test statistic) do we
                % need to also exclude if tract is connected to too many VTA
                gval{side}(Nmap>((1-(obj.statsettings.connthreshold/100))*length(gpatsel)),gpatsel)=nan;
            end
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

if strcmp(obj.threshstrategy,'Unthresholded')
    if strcmp(obj.drawTool,'fiberfiltering')
        for side=1:numel(gval)
            usedidx{group,side} = find(isfinite(vals{group,side}));
            fibcell{group,side} = obj.results.fiberfiltering.(ea_conn2connid(obj.calcsettings.fibfilt_connectome)).fibcell{side}(usedidx{group,side});
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
                fibcell{group,side} = obj.results.fiberfiltering.(ea_conn2connid(obj.calcsettings.fibfilt_connectome)).fibcell{side}(usedidx{group,side});
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
