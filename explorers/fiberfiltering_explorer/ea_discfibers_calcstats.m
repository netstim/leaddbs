function [vals,fibcell,usedidx] = ea_discfibers_calcstats(obj,patsel,Iperm)

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

% fiber values can be sigmoid transform
switch obj.statsettings.stimulationmodel
    case 'Sigmoid Field'
        if obj.connectivity_type == 2
            fibsval = obj.results.(ea_conn2connid(obj.connectome)).('PAM_probA').fibsval;
        else
            fibsval_raw = obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval;
            fibsval = fibsval_raw;  % initialize
            for side = 1:size(fibsval_raw,2)
                fibsval{1,side}(:,:) = ea_SigmoidFromEfield(fibsval_raw{1,side}(:,:));
            end
        end
    otherwise
        fibsval = cellfun(@full, obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval, 'Uni', 0);
end

if size(I,2)==1 % 1 entry per patient, not per electrode
    I=[I,I]; % both sides the same;
end

if obj.mirrorsides
    I=[I;I];
end

%% Centrally deal with Covariates here by regressing them out from target variable:
if ~isempty(obj.covars)
    for i=1:length(obj.covars)
        if obj.mirrorsides
            covars{i} = [obj.covars{i};obj.covars{i}];
        else
            covars{i} = obj.covars{i};
        end
        % directly regress out covars from variable here:
        for side=1:2
            if size(covars{i},2)==1 % single entry, use the same for both sides
                I(:,side)=ea_resid(covars{i},I(:,side));
            elseif size(covars{i},2)==2 % two entries in covars, use each for the respective side:
                I(:,side)=ea_resid(covars{i}(:,side),I(:,side));
            end
        end
    end
end

if strcmp(obj.multitractmode,'Split & Color By PCA')
    % prep PCA: here, we need to get scores for all the patients. But PCA should
    % remain based on the patients selected for analysis 
    
    % variables for the selected patients 
    for i=1:length(obj.subscore.vars)
        selected_subscores{i} = obj.subscore.vars{i}(obj.patientselection);
    end
    selected_subscores = cell2mat(selected_subscores);
    subvars=ea_nanzscore(selected_subscores);
    
    % PCA - get PC scores for selected patients
    [coeff,score,latent,tsquared,explained,mu]=pca(subvars,'rows','complete');
    % [coeff,score,latent,tsquared,explained,mu]=pca(subvars,'rows','pairwise'); %pca
    
    if isempty(score) || ea_isnan(score,'any')
        score=nan(length(obj.responsevar),obj.numpcs);
    end

    obj.subscore.pcavars=cell(obj.numpcs,1);

    for pc=1:obj.numpcs
        obj.subscore.pcavars{pc}(obj.patientselection,1)=score(:,pc); %pca variables -> pca components, location of first subscore is replaced by first pc
    end
    if ~isfield(obj.subscore,'pcacolors')
        obj.subscore.pcacolors=ea_color_wes('all'); % assign some random colors.
    end

    % now use PCA weights to get PC scores for the non-selected patients 
    patientnonsel = logical(ones(1, length(obj.allpatients)));
    patientnonsel(obj.patientselection) = 0; 

    for i=1:length(obj.subscore.vars)
        nonsel_subscores{i} = obj.subscore.vars{i}(patientnonsel);
    end
    nonsel_subscores = cell2mat(nonsel_subscores);
    % pseudo zscore - use mean and sd of selected patients to keep same "scale"
    for ci = 1:size(selected_subscores, 2)
        datawonan = selected_subscores(:, ci);
        datawonan = datawonan(~isnan(datawonan));
        datamean(ci) = mean(datawonan);
        datasd(ci) = std(datawonan);
    end
    nonsel_subvars = ( nonsel_subscores - repmat(datamean, size(nonsel_subscores, 1), 1) ) ...
        ./ repmat(datasd, size(nonsel_subscores, 1), 1);
    
    % multiply clinical scores by weights
    for pc=1:obj.numpcs
        obj.subscore.pcavars{pc}(patientnonsel,1)= nonsel_subvars*coeff(:,pc);
    end

    % save the PCA coefficients for later 
    obj.subscore.pcacoeff = coeff; 

end

switch obj.multitractmode
    case 'Split & Color By Group'
        groups = unique(obj.M.patient.group(obj.patientselection))';
        dogroups = 1;
        dosubscores = 0;
    case 'Split & Color By Subscore'
        if ~isempty(obj.subscore.vars) %this will be empty when user
            %initializes the split by subscore button
           groups = 1:length(obj.subscore.vars);
           dosubscores = 1;
           dogroups = 0;
        else
            obj.multitractmode = 'Single Tract Analysis';
            groups = 1;
            dogroups = 0;
            dosubscores = 0;
        end
    case 'Split & Color By PCA'
        groups = 1:length(obj.subscore.pcavars);
        dosubscores = 1;
        dogroups = 0;
    otherwise
        groups=1;
        dogroups = 0;
        dosubscores = 0;
end

for group=groups
    gfibsval=fibsval; %refresh fibsval
    if dogroups
        groupspt=find((obj.M.patient.group)'==group);
        gpatsel=groupspt(ismember(groupspt,patsel));
    elseif dosubscores
        gpatsel=patsel;
    else
        gpatsel=patsel;
    end
    if obj.mirrorsides
        gpatsel=[gpatsel,gpatsel + length(obj.allpatients)];
    end
    if dosubscores
        switch obj.multitractmode
            case 'Split & Color By Subscore'
                I = obj.subscore.vars{group};
            case 'Split & Color By PCA'
                if ~exist('Iperm','var')
                    I = obj.subscore.pcavars{group};
                else 
                    % recompute permuted PC scores based on permuted
                    % clinical vars - keep right dimensions for I but use
                    % only selected patients for zscores - patientsel
                    % applied to I later
                    I = nan(size(Iperm)); 
                    I(obj.patientselection,:) = ea_nanzscore(Iperm(obj.patientselection,:))*coeff; 
                    I = I(:, group); 
                end    
        end
        if size(I,2)==1 % 1 entry per patient, not per electrode
            I=[I,I]; % both sides the same;
        end
        if obj.mirrorsides
            I=[I;I];
        end
    end

    for side=1:numel(gfibsval)

        %% Step 1: Prefiltering. This part does not use improvements but simply selects potential tracts that fullfil criteria (conn sliders):

        % check connthreshold
        if obj.runwhite || strcmp(obj.statsettings.stattest,'N-Map')
            Nmap=sum(gfibsval{side}(:,gpatsel),2);
        else
            switch obj.statsettings.stimulationmodel
                case 'VTA'
                    Nmap=sum(gfibsval{side}(:,gpatsel),2);
%                 case 'Sigmoid Field'
%                     if strcmp(ea_method2methodid(obj), 'spearman_5peak') || strcmp(ea_method2methodid(obj), 'spearman_peak')
%                         % 0.5 V / mm -> 0.5 probability
%                         Nmap=sum((gfibsval{side}(:,gpatsel)>obj.statsettings.efieldthreshold/1000.0),2);
%                     else
%                         Nmap=sum((gfibsval{side}(:,gpatsel)>obj.statsettings.efieldthreshold),2);
%                     end
                case 'Sigmoid Field'
                    Nmap=sum((gfibsval{side}(:,gpatsel)>obj.statsettings.efieldthreshold),2);
                case 'Electric Field'
                    Nmap=sum((gfibsval{side}(:,gpatsel)>obj.statsettings.efieldthreshold),2);
            end
        end
        % remove fibers that are not connected to enough VTAs/Efields or connected
        % to too many VTAs (connthreshold slider)
        if ~obj.runwhite
            gfibsval{side}(Nmap<((obj.statsettings.connthreshold/100)*length(gpatsel)),gpatsel)=nan;
            if strcmp(obj.statsettings.stimulationmodel,'VTA')
                % only in case of VTAs (given two-sample-t-test statistic) do we
                % need to also exclude if tract is connected to too many VTAs:
                gfibsval{side}(Nmap>((1-(obj.statsettings.connthreshold/100))*length(gpatsel)),gpatsel)=nan;
            end
        end

        % init outputvars
        vals{group,side}=nan(size(gfibsval{side},1),1);
        if obj.showsignificantonly
            pvals{group,side}=vals{group,side};
        end
        if obj.runwhite
            vals{group,side} = Nmap/length(gpatsel);
            vals{group,side}(vals{group,side}==0)=nan;
        else
            nonempty=sum(gfibsval{side}(:,gpatsel),2,'omitnan')>0;
            nonemptyidx=find(nonempty);

            valsin=gfibsval{side}(nonempty,gpatsel);
            outcomein=I(gpatsel,side);

            if strcmp(obj.statsettings.statfamily, 'Correlations')
                disp(['Calculating ' obj.statsettings.stattest ' correlation for side ' num2str(side) '...']);
            else
                disp(['Calculating ' obj.statsettings.stattest ' for side ' num2str(side) '...']);
            end

            stattests=ea_explorer_statlist(obj.responsevar);

            [is,idx]=ismember(obj.statsettings.stattest,stattests.name);
            if ~is
                ea_error(['Function for test ',obj.statsettings.stattest,' missing.']);
            end

            %% Step 2: Fiberfiltering. This part filters fibers based on outcome variable (except for descriptive tests):

            %this following line calls the actual statistical test:
            if ~isempty(valsin) && ~isempty(outcomein)
                [valsout,psout]=feval(stattests.file(idx),valsin,outcomein,obj.statsettings.H0); % apply test
                vals{group,side}(nonemptyidx)=valsout;
                if exist('pvals','var')
                    pvals{group,side}(nonemptyidx)=psout;
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
end

% close group loop to test for significance across all tests run:
if ~obj.runwhite
    if obj.showsignificantonly
        vals=ea_corrsignan(vals,pvals,obj);
    end
end

% Clean up non-finite values from fibcell and vals
for group=groups
    for side=1:numel(gfibsval)
        usedidx{group,side} = find(isfinite(vals{group,side}));
        fibcell{group,side} = obj.results.(ea_conn2connid(obj.connectome)).fibcell{side}(usedidx{group,side});
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

unthresholdedVals = vals; % Need to keep this original vals when calculating the (same) threshold to both sides. 

% Thresholding
for group=groups
    for side=1:numel(gfibsval)
        switch obj.threshstrategy
            case 'Fixed Amount' % here we want to create thresholds for each side separately.
                posvals = sort(vals{group,side}(vals{group,side}>0),'descend');
                negvals = sort(vals{group,side}(vals{group,side}<0),'ascend');
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
            vals{group,side}(remove)=[];
            fibcell{group,side}(remove)=[];
            usedidx{group,side}(remove)=[];
        end
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
