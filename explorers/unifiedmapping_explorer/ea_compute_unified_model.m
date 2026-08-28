function [Ihat,Ihat_train_global,val_struct] = ea_compute_unified_model(numTestIt, obj, fibsval, Ihat, Ihat_train_global, patientsel, training, test, Iperm,vals,usedidx)

% this needs to be edited when external model is brought back
% if obj.useExternalModel == true
%         S = load(obj.ExternalModelFile);
%         if ~strcmp(S.connectome,ea_conn2connid(obj.connectome))
%             waitfor(msgbox('The chosen fibfilt model was computed for another connectome! See terminal'));
%             disp('Model for connectome: ')
%             disp(S.connectome)
%             return
%         end
% 
%         if obj.calcsettings.connectivity_type ~= S.conn_type
%             waitfor(msgbox('The connectivity methods of imported and current model are different! See terminal'));
%             disp('Connectivity type of imported model: ')
%             disp(obj.calcsettings.connectivity_type)
%             return     
%         end
% 
%         vals_connected = cell(size(S.vals_all,1),size(S.vals_all,2));
%         for voter = 1:size(vals_connected,1)
%             for side=1:size(vals_connected,2)
%                 try
%                     switch obj.calcsettings.connectivity_type
%                         case 2
%                             if size(S.vals_all,2) == 2
%                                 vals_connected{voter,side} = S.vals_all{voter,side}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{side});
%                             else
%                                 vals_connected{voter,side} = S.vals_all{voter,1}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{side});
%                             end
%                         otherwise
%                             if size(S.vals_all,2) == 2
%                                 vals_connected{voter,side} = S.vals_all{voter,side}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{side});
%                             else
%                                 vals_connected{voter,side} = S.vals_all{voter,1}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{side});
%                             end
%                     end
%                 catch
%                     ea_warndlg("Connectivity indices were not stored. Please recalculate.");
%                     return
%                 end
%             end
%         end
%     end

    % Sigmoid fiber values are calculated and stored with the connectivity
    % results. Do not transform conventional E-field fiber values here.
    if strcmp(obj.statsettings.stimulationmodel, 'Sigmoid Field')
        if obj.calcsettings.connectivity_type == 2
            fibsval = obj.results.fiberfiltering.(ea_conn2connid(obj.calcsettings.fibfilt_connectome)).('PAM_probA').fibsval;
        end
    end

    if ~exist('Iperm', 'var') || isempty(Iperm)
        
        [vals,fibcell,usedidx] = ea_unifiedmapping_calcstats(obj, patientsel(training));
    else
        [vals,fibcell,usedidx] = ea_unifiedmapping_calcstats(obj, patientsel(training), Iperm);
    end

    if obj.cvlivevisualize
        obj.draw(vals,fibcell,usedidx)
        drawnow;
    end
    % if obj.useExternalModel == true
    %    [vals,fibcell,usedidx] = ea_discfibers_loadModel_calcstats(obj, vals_connected);
    % else
    %    [vals,fibcell,usedidx] = ea_unified_discfibers_calcstats(obj, patientsel(training));
    %end
    % if no fibers were selected for the permutation iteration,
    % assign dummies that will have r = 0
    if isempty(vals{1}) && isempty(vals{2}) && exist('Iperm', 'var') && ~isempty(Iperm)
        for voter=1:size(vals,1)
            for side=1:size(vals,2)
                Ihat(test,side,voter) = 42;
                Ihat_train_global(numTestIt,training,side,voter) = 42;
            end
        end

        val_struct.vals=vals;
        val_struct.usedidx=usedidx;
        try    
            val_struct.fibcell=fibcell; 
        end % fibcell not always supplied.

        return
    end

    if size(obj.responsevar,2)>1
        lateral_score = true;
    else
        lateral_score = false;
    end
    if strcmp(obj.drawTool,'networkmapping')
        lateral_score = false;
    end
    switch obj.modelNormalization
        case 'z-score'
            for s=1:length(vals)
                vals{s}=ea_nanzscore(vals{s});
            end
        case 'van Albada 2007'
            for s=1:length(vals)
                vals{s}=ea_normal(vals{s});
            end
    end

    for voter=1:size(vals,1)
        % retrieve indices of all connected fibers and pivotal fibers
        if lateral_score == false && size(vals,2) == 2
            vals_flat = vertcat(vals{voter,:});
            % just concatenate values from both hemispheres
        else
            vals_flat=vals{voter};
        end
            orig_model = cell(1,2);
            switch obj.drawTool
                case 'sweetspotmapping'
                    orig_model{1} = obj.results.sweetspotmapping.efield{1}';
                    if size(obj.results.sweetspotmapping.efield,1)==2
                        orig_model{2} = obj.results.sweetspotmapping.efield{2}';
                    end
                    % Apply the training sweet spot's missing-value threshold
                    % to raw E-fields before any sigmoid transformation.
                    if ismember(obj.statsettings.stimulationmodel, ...
                            {'Electric Field', 'Sigmoid Field'}) && ...
                            isfield(obj.statsettings, 'nanthreshold') && ...
                            ~isempty(obj.statsettings.nanthreshold)
                        for modelSide = 1:numel(orig_model)
                            if ~isempty(orig_model{modelSide})
                                orig_model{modelSide}( ...
                                    orig_model{modelSide} <= ...
                                    obj.statsettings.nanthreshold) = nan;
                            end
                        end
                    end

                    % Training statistics use sigmoid-transformed E-fields,
                    % so prediction must use the same representation.
                    if strcmp(obj.statsettings.stimulationmodel, 'Sigmoid Field')
                        for modelSide = 1:numel(orig_model)
                            if ~isempty(orig_model{modelSide})
                                orig_model{modelSide} = ...
                                    ea_SigmoidFromEfield(orig_model{modelSide});
                            end
                        end
                    end
                    orig_model_flat = vertcat(orig_model{:});
                case 'fiberfiltering'
                    orig_model{1} = fibsval{1,1}(usedidx{voter,1},patientsel);
                    if size(vals,2) == 2
                        orig_model{2} = fibsval{1,2}(usedidx{voter,2},patientsel);
                    end
                    orig_model_flat = vertcat(orig_model{:});
                case 'networkmapping'
                    orig_model_flat = (full(obj.results.networkmapping.(ea_conn2connid(obj.calcsettings.netmap_connectome)).connval))';
            end
  
        for side=1:size(vals,2)
            % if not lateral, compute Ihat for both hemispheres at the same time
            if ~isempty(vals{voter,side})
                % Sweet spots are scored per hemisphere and combined later
                % by the caller. All other mapping modes retain their
                % original lateral/non-lateral behavior.
                if strcmp(obj.drawTool, 'sweetspotmapping')
                    side_model = orig_model{side}(:,patientsel);
                    Ihat(test,side,voter) = ...
                        sweetspotSidePrediction( ...
                            vals{voter,side}, side_model(:,test), obj);
                    Ihat_train_global(numTestIt,training,side,voter) = ...
                        sweetspotSidePrediction( ...
                            vals{voter,side}, side_model(:,training), obj);
                    continue
                elseif strcmp(obj.drawTool, 'fiberfiltering')
                    combine_sides_before_prediction = ...
                        lateral_score == false;
                    side_model = fibsval{1,side}( ...
                        usedidx{voter,side},patientsel);
                else
                    combine_sides_before_prediction = ...
                        lateral_score == false;
                    side_model = [];
                end

                switch obj.statsettings.stimulationmodel % also differentiate between methods in the prediction part.
                    case {'VTA'} % VTAs
                        switch lower(obj.basepredictionon)
                            case 'mean of scores'
                                if combine_sides_before_prediction
                                    Ihat_all = ea_nanmean(vals_flat.*orig_model_flat,1);
                                    Ihat(test,1, voter) = Ihat_all(test);
                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(orig_model_flat(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,1, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA
                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);
                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_nanmean(vals{voter,side}.*side_model(:,test),1);
                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(side_model(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,side, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA
                                    Ihat_train_global(numTestIt,training,side,voter) = ea_nanmean(vals{voter,side}.*side_model(:,training),1);
                                end

                            case 'sum of scores'
                                if combine_sides_before_prediction
                                    Ihat_all = ea_nansum(vals_flat.*orig_model_flat,1);
                                    Ihat(test,1, voter) = Ihat_all(test);
                                    
                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(orig_model_flat(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,1, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA

                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);
                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_nansum(vals{voter,side}.*side_model(:,test),1);
                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(side_model(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,side, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA
                                    Ihat_train_global(numTestIt,training,side,voter) = ea_nansum(vals{voter,side}.*side_model(:,training),1);
                                 
                                end

                            case 'peak of scores'
                                if combine_sides_before_prediction
                                    Ihat_all = ea_discfibers_getpeak(vals_flat.*orig_model_flat, obj.posvisible, obj.negvisible, 'peak');
                                    Ihat(test,1, voter) = Ihat_all(test);
                                    
                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(orig_model_flat(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,1, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA

                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);


                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_discfibers_getpeak(vals{voter,side}.*side_model(:,test), obj.posvisible, obj.negvisible, 'peak');
                                    
                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(side_model(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,side, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA

                                    Ihat_train_global(numTestIt,training,side,voter) = ea_discfibers_getpeak(vals{voter,side}.*side_model(:,training), obj.posvisible, obj.negvisible, 'peak');
                                   
                                end
                            case 'peak 5% of scores'
                                if combine_sides_before_prediction
                                    Ihat_all = ea_discfibers_getpeak(vals_flat.*orig_model_flat, obj.posvisible, obj.negvisible, 'peak5');
                                    Ihat(test,1, voter) = Ihat_all(test);

                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(orig_model_flat(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,1, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA

                                    Ihat(test,2, voter) = Ihat(test,1, voter);

                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);

                                    break % both sides are already filled out!
                                else
                                    if strcmp(obj.drawTool, 'sweetspotmapping')
                                        ihatvals=vals{voter,side}.*side_model;
                                    else
                                        % Preserve the existing VTA fiber-filtering
                                        % behavior for multi-voter models.
                                        ihatvals=vals{1,side}.*side_model;
                                    end
                                    Ihat(test,side,voter) = ea_discfibers_getpeak(ihatvals(test), obj.posvisible, obj.negvisible, 'peak5');
                                    
                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(side_model(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,side, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA
                                    
                                    Ihat_train_global(numTestIt,training,side,voter) = ea_discfibers_getpeak(ihatvals(training), obj.posvisible, obj.negvisible, 'peak5');
                                end
                        end
                    case {'Electric Field', 'Sigmoid Field'} % efields
                        switch lower(obj.basepredictionon)
                            case 'profile of scores: spearman'
                                if combine_sides_before_prediction
                                    Ihat_all = corr(vals_flat,orig_model_flat,'rows','pairwise','type','spearman');
                                    Ihat(test,1, voter) = Ihat_all(test);

                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(orig_model_flat(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,1, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA

                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_all(training);

                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = corr(vals{voter,side},side_model(:,test),'rows','pairwise','type','spearman');
                                    
                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(side_model(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,side, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA
                                    Ihat_train_global(numTestIt,training,side,voter) = corr(vals{voter,side},side_model(:,training),'rows','pairwise','type','spearman');
                                end
                            case 'profile of scores: pearson'
                                if combine_sides_before_prediction
                                    Ihat_all = corr(vals_flat,orig_model_flat,'rows','pairwise','type','pearson');
                                    Ihat(test,1, voter) = Ihat_all(test);

                                    testidx=find(test);
                                    allzerotestidx = testidx(sum(orig_model_flat(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,1, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA


                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);

                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = corr(vals{voter,side},side_model(:,test),'rows','pairwise','type','pearson');
                                    
                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(side_model(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,side, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA

                                    Ihat_train_global(numTestIt,training,side,voter) = corr(vals{voter,side},side_model(:,training),'rows','pairwise','type','pearson');
                                end
                            case 'profile of scores: bend'
                                if combine_sides_before_prediction
                                    Ihat_all = ea_bendcorr(vals_flat,orig_model_flat);
                                    Ihat(test,1, voter) = Ihat_all(test);

                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(orig_model_flat(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,1, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA

                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);

                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_bendcorr(vals{voter,side},side_model(:,test));
                                    
                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(side_model(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,side, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA
                                    
                                    Ihat_train_global(numTestIt,training,side,voter) = ea_bendcorr(vals{voter,side},side_model(:,training));
                                    
                                end
                            case 'mean of scores'
                                if combine_sides_before_prediction
                                    Ihat_all = ea_nanmean(vals_flat.*orig_model_flat,1);
                                    Ihat(test,1, voter) = Ihat_all(test);
                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(orig_model_flat(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,1, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA
                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_all(training);

                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_nanmean(vals{voter,side}.*side_model(:,test),1);
                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(side_model(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,side, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA
                                    Ihat_train_global(numTestIt,training,side,voter) = ea_nanmean(vals{voter,side}.*side_model(:,training),1);
                                end
                            case 'sum of scores'
                                if combine_sides_before_prediction
                                    Ihat_all = ea_nansum(vals_flat.*orig_model_flat,1);
                                    Ihat(test,1, voter) = Ihat_all(test);

                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(orig_model_flat(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,1, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA

                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);

                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_nansum(vals{voter,side}.*side_model(:,test),1);
                                    
                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(side_model(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,side, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA

                                    Ihat_train_global(numTestIt,training,side,voter) = ea_nansum(vals{voter,side}.*side_model(:,training),1);
                                   
                                end
                            case 'peak of scores'
                                if combine_sides_before_prediction
                                    Ihat_all = ea_discfibers_getpeak(vals_flat.*orig_model_flat, obj.posvisible, obj.negvisible, 'peak');
                                    Ihat(test,1, voter) = Ihat_all(test);

                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(orig_model_flat(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,1, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA

                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);

                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_discfibers_getpeak(vals{voter,side}.*side_model(:,test), obj.posvisible, obj.negvisible, 'peak');
                                   
                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(side_model(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,side, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA

                                    Ihat_train_global(numTestIt,training,side,voter) = ea_discfibers_getpeak(vals{voter,side}.*side_model(:,training), obj.posvisible, obj.negvisible, 'peak');
                                end
                            case 'peak 5% of scores'
                                if combine_sides_before_prediction
                                    Ihat_all = ea_discfibers_getpeak(vals_flat.*orig_model_flat, obj.posvisible, obj.negvisible, 'peak5');
                                    Ihat(test,1, voter) = Ihat_all(test);

                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(orig_model_flat(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,1, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA

                                    Ihat(test,2, voter) = Ihat(test,1, voter);

                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);

                                    break % both sides are already filled out!
                                else
                                    ihatvals=vals{voter,side}.*side_model;
                                    Ihat(test,side,voter) = ea_discfibers_getpeak(ihatvals(test), obj.posvisible, obj.negvisible, 'peak5');
                                    
                                    testidx=find(test);
                                    allzerotestidx=testidx(sum(side_model(:,test),1,'omitnan') == 0);
                                    Ihat(allzerotestidx,side, voter) = nan; % set Ihats to nan if there is no overlap with even a single VTA

                                    Ihat_train_global(numTestIt,training,side,voter) = ea_discfibers_getpeak(ihatvals(training), obj.posvisible, obj.negvisible, 'peak5');
                                end
                        end
                end
            end
        end

    end
    val_struct.vals=vals;
    val_struct.usedidx=usedidx;
try    
    val_struct.fibcell=fibcell; end % fibcell not always supplied.

end

function prediction = sweetspotSidePrediction(vals, efield, obj)
% Score one sweet-spot hemisphere without fiber-specific overlap handling.

switch lower(obj.basepredictionon)
    case 'profile of scores: spearman'
        prediction = corr(vals, efield, ...
            'rows', 'pairwise', 'type', 'spearman');
    case 'profile of scores: pearson'
        prediction = corr(vals, efield, ...
            'rows', 'pairwise', 'type', 'pearson');
    case 'profile of scores: bend'
        prediction = ea_bendcorr(vals, efield);
    case 'mean of scores'
        prediction = ea_nanmean(vals .* efield, 1);
    case 'sum of scores'
        prediction = ea_nansum(vals .* efield, 1);
    case 'peak of scores'
        prediction = ea_discfibers_getpeak( ...
            vals .* efield, obj.posvisible, obj.negvisible, 'peak');
    case 'peak 5% of scores'
        prediction = ea_discfibers_getpeak( ...
            vals .* efield, obj.posvisible, obj.negvisible, 'peak5');
    otherwise
        error('Unknown prediction method: %s', obj.basepredictionon);
end
end
