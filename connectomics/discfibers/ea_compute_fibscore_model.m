function [Ihat,Ihat_train_global,val_struct,actualimprovs] = ea_compute_fibscore_model(numTestIt,adj_scaler, obj, fibsval, Ihat, Ihat_train_global, patientsel, training, test, Iperm)

    if obj.useExternalModel == true
        S = load(obj.ExternalModelFile);
        if ~strcmp(S.connectome,ea_conn2connid(obj.connectome))
            waitfor(msgbox('The chosen fibfilt model was computed for another connectome! See terminal'));
            disp('Model for connectome: ')
            disp(S.connectome)
            return
        end

        if obj.connectivity_type ~= S.conn_type
            waitfor(msgbox('The connectivity methods of imported and current model are different! See terminal'));
            disp('Connectivity type of imported model: ')
            disp(obj.connectivity_type)
            return     
        end

        vals_connected = cell(size(S.vals_all,1),size(S.vals_all,2));
        for voter = 1:size(vals_connected,1)
            for side=1:size(vals_connected,2)
                try
                    switch obj.connectivity_type
                        case 2
                            if size(S.vals_all,2) == 2
                                vals_connected{voter,side} = S.vals_all{voter,side}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{side});
                            else
                                vals_connected{voter,side} = S.vals_all{voter,1}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{side});
                            end
                        otherwise
                            if size(S.vals_all,2) == 2
                                vals_connected{voter,side} = S.vals_all{voter,side}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{side});
                            else
                                vals_connected{voter,side} = S.vals_all{voter,1}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{side});
                            end
                    end
                catch
                    ea_warndlg("Connectivity indices were not stored. Please recalculate.");
                    return
                end
            end
        end
    end

    if ~exist('Iperm', 'var') || isempty(Iperm)
        if obj.cvlivevisualize
            if obj.useExternalModel == true
              [vals,fibcell,usedidx] = ea_discfibers_loadModel_calcstats(obj, vals_connected);
            else
              [vals,fibcell,usedidx] = ea_discfibers_calcstats(obj, patientsel(training));
            end

            %[vals,fibcell,usedidx] = ea_discfibers_calcstats(obj, patientsel(training));
            obj.draw(vals,fibcell,usedidx)
            %obj.draw(vals,fibcell);
            drawnow;
        else
            if obj.useExternalModel == true
               [vals,~,usedidx] = ea_discfibers_loadModel_calcstats(obj, vals_connected);
            else
               [vals,~,usedidx] = ea_discfibers_calcstats(obj, patientsel(training));
            end
            %[vals,~,usedidx] = ea_discfibers_calcstats(obj, patientsel(training));
        end
    else
        if obj.cvlivevisualize
            [vals,fibcell,usedidx] = ea_discfibers_calcstats(obj, patientsel(training), Iperm);
            obj.draw(vals,fibcell,usedidx)
            %obj.draw(vals,fibcell);
            drawnow;
        else
            [vals,~,usedidx] = ea_discfibers_calcstats(obj, patientsel(training), Iperm);
        end
    end

    % if no fibers were selected for the permutation iteration, 
    % assign dummies that will have r = 0
    if isempty(vals{1}) && isempty(vals{2}) && exist('Iperm', 'var') && ~isempty(Iperm)
        for voter=1:size(vals,1)
            for side=1:size(vals,2)
                actualimprovs{voter,side} = 42;
                Ihat(test,side,voter) = 42;
                Ihat_train_global(numTestIt,training,side,voter) = 42;
            end
        end

        val_struct.vals=vals;
        val_struct.usedidx=usedidx;
        try    val_struct.fibcell=fibcell; end % fibcell not always supplied.

        return
    end

    if size(obj.responsevar,2)>1
        lateral_score = true;
    else
        lateral_score = false;
    end
    if obj.CleartuneOptim
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
        if isstruct(obj.ADJ)
            if obj.connectivity_type == 2
                global_val_ind = obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{1}(usedidx{voter,1});
                global_conn_ind = obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{1}(1:end);
                global_val_ind_lh = obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{2}(usedidx{voter,2});
                global_conn_ind_lh = obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{2}(1:end);
            else
                global_val_ind = obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{1}(usedidx{voter,1});
                global_conn_ind = obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{1}(1:end);
                global_val_ind_lh = obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{2}(usedidx{voter,2});
                global_conn_ind_lh = obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{2}(1:end);

            end
        end

        for side=1:size(vals,2)

            % if not lateral, compute Ihat for both hemispheres at the same time
            if lateral_score == false
                % just concatenate values from both hemispheres
                fibsval_usedidx = cell(1,2);
                fibsval_usedidx{1} = fibsval{1,1}(usedidx{voter,1},patientsel);
                if size(vals,2)>1
                    fibsval_usedidx{2} = fibsval{1,2}(usedidx{voter,2},patientsel);
                end
                vals_flat = vertcat(vals{voter,:});
                fibsval_usedidx_flat = vertcat(fibsval_usedidx{:});
            end

            if ~isempty(vals{voter,side})
                switch obj.statmetric % also differentiate between methods in the prediction part.
                    case {'One-Sample Tests / VTAs / PAM (OSS-DBS)','Two-Sample T-Tests / VTAs (Baldermann 2019) / PAM (OSS-DBS)', 'Proportion Test (Chi-Square) / VTAs (binary vars)'...
                            'Binomial Tests / VTAs (binary vars)', 'Plain Connections'} % VTAs
                        switch lower(obj.basepredictionon)
                            case 'mean of scores'
                                if lateral_score == false
                                    Ihat_all = ea_nanmean(vals_flat.*fibsval_usedidx_flat,1);
                                    Ihat(test,1, voter) = Ihat_all(test);
                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);

                                    if isstruct(obj.ADJ)
                                        % compute fiberscores for the rest of the
                                        % connected fibers (weighted by the correposponding vals)
                                        Ihat_ADJ_rh = vals{voter,1}.*(obj.ADJ.ADJ(global_val_ind,global_conn_ind)*fibsval{1,1}(1:end,patientsel));
                                        Ihat_ADJ_lh = vals{voter,2}.*(obj.ADJ.ADJ(global_val_ind_lh,global_conn_ind_lh)*fibsval{1,2}(1:end,patientsel));
                                        Ihat_ADJ_comb = [Ihat_ADJ_rh;Ihat_ADJ_lh];
                                        % compute the score for both
                                        % hemispheres together, adj_scaler defines the extent of 'smoothing'
                                        Ihat_ADJ_test = transpose(adj_scaler*ea_nanmean(Ihat_ADJ_comb(1:end,test),1));
                                        Ihat_ADJ_training = adj_scaler*ea_nanmean(Ihat_ADJ_comb(1:end,training),1);
                                        Ihat(test,1,voter) = Ihat(test,1,voter) + Ihat_ADJ_test;
                                        Ihat(test,2, voter) = Ihat(test,1, voter);
                                        Ihat_train_global(numTestIt,training,1, voter) = Ihat_train_global(numTestIt,training,1, voter) + Ihat_ADJ_training;
                                        Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1, voter);
                                    end

                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                    Ihat_train_global(numTestIt,training,side,voter) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                    if isstruct(obj.ADJ)
                                        disp('Adjacency matrix for lateral symptoms is currently not supported')
                                    end
                                end

                            case 'sum of scores'
                                if lateral_score == false
                                    Ihat_all = ea_nansum(vals_flat.*fibsval_usedidx_flat,1);
                                    Ihat(test,1, voter) = Ihat_all(test);
                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);

                                    if isstruct(obj.ADJ)
                                        % compute fiberscores for the rest of the
                                        % connected fibers (weighted by the correposponding vals)
                                        Ihat_ADJ_rh = vals{voter,1}.*(obj.ADJ.ADJ(global_val_ind,global_conn_ind)*fibsval{1,1}(1:end,patientsel));
                                        Ihat_ADJ_lh = vals{voter,2}.*(obj.ADJ.ADJ(global_val_ind_lh,global_conn_ind_lh)*fibsval{1,2}(1:end,patientsel));
                                        Ihat_ADJ_comb = [Ihat_ADJ_rh;Ihat_ADJ_lh];
                                        Ihat_ADJ_test = transpose(adj_scaler*ea_nansum(Ihat_ADJ_comb(1:end,test),1));
                                        Ihat_ADJ_training = adj_scaler*ea_nansum(Ihat_ADJ_comb(1:end,training),1);
                                        Ihat(test,1,voter) = Ihat(test,1,voter) + Ihat_ADJ_test;
                                        Ihat(test,2, voter) = Ihat(test,1, voter);
                                        Ihat_train_global(numTestIt,training,1, voter) = Ihat_train_global(numTestIt,training,1, voter) + Ihat_ADJ_training;
                                        Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1, voter);
                                    end

                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                    Ihat_train_global(numTestIt,training,side,voter) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                    if isstruct(obj.ADJ)
                                        disp('Adjacency matrix for lateral symptoms is currently not supported')
                                    end
                                end

                            case 'peak of scores'
                                if lateral_score == false
                                    Ihat_all = ea_discfibers_getpeak(vals_flat.*fibsval_usedidx_flat, obj.posvisible, obj.negvisible, 'peak');
                                    Ihat(test,1, voter) = Ihat_all(test);
                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);

                                    if isstruct(obj.ADJ)
                                        % compute fiberscores for the rest of the
                                        % connected fibers (weighted by the correposponding vals)
                                        Ihat_ADJ_rh = vals{voter,1}.*(obj.ADJ.ADJ(global_val_ind,global_conn_ind)*fibsval{1,1}(1:end,patientsel));
                                        Ihat_ADJ_lh = vals{voter,2}.*(obj.ADJ.ADJ(global_val_ind_lh,global_conn_ind_lh)*fibsval{1,2}(1:end,patientsel));
                                        Ihat_ADJ_comb = [Ihat_ADJ_rh;Ihat_ADJ_lh];
                                        Ihat_ADJ_test = transpose(adj_scaler*ea_discfibers_getpeak(Ihat_ADJ_comb(1:end,test), obj.posvisible, obj.negvisible, 'peak'));
                                        Ihat_ADJ_training = adj_scaler*ea_discfibers_getpeak(Ihat_ADJ_comb(1:end,training), obj.posvisible, obj.negvisible, 'peak');
                                        Ihat(test,1,voter) = Ihat(test,1,voter) + Ihat_ADJ_test;
                                        Ihat(test,2, voter) = Ihat(test,1, voter);
                                        Ihat_train_global(numTestIt,training,1, voter) = Ihat_train_global(numTestIt,training,1, voter) + Ihat_ADJ_training;
                                        Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1, voter);
                                    end

                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_discfibers_getpeak(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)), obj.posvisible, obj.negvisible, 'peak');
                                    Ihat_train_global(numTestIt,training,side,voter) = ea_discfibers_getpeak(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)), obj.posvisible, obj.negvisible, 'peak');
                                    if isstruct(obj.ADJ)
                                        disp('Adjacency matrix for lateral symptoms is currently not supported')
                                    end
                                end
                            case 'peak 5% of scores'
                                if lateral_score == false
                                    Ihat_all = ea_discfibers_getpeak(vals_flat.*fibsval_usedidx_flat, obj.posvisible, obj.negvisible, 'peak5');
                                    Ihat(test,1, voter) = Ihat_all(test);
                                    Ihat(test,2, voter) = Ihat(test,1, voter);

                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);

                                    if isstruct(obj.ADJ)
                                        Ihat_ADJ_rh = vals{voter,1}.*(obj.ADJ.ADJ(global_val_ind,global_conn_ind)*fibsval{1,1}(1:end,patientsel));
                                        Ihat_ADJ_lh = vals{voter,2}.*(obj.ADJ.ADJ(global_val_ind_lh,global_conn_ind_lh)*fibsval{1,2}(1:end,patientsel));
                                        Ihat_ADJ_comb = [Ihat_ADJ_rh;Ihat_ADJ_lh];

                                        Ihat_ADJ_test = transpose(adj_scaler*ea_discfibers_getpeak(Ihat_ADJ_comb(test), obj.posvisible, obj.negvisible, 'peak5'));
                                        Ihat_ADJ_training = adj_scaler*ea_discfibers_getpeak(Ihat_ADJ_comb(training), obj.posvisible, obj.negvisible, 'peak5');
                                        Ihat(test,1,voter) = Ihat(test,1,voter) + Ihat_ADJ_test;
                                        Ihat(test,2, voter) = Ihat(test,1, voter);
                                        Ihat_train_global(numTestIt,training,1, voter) = Ihat_train_global(numTestIt,training,1, voter) + Ihat_ADJ_training;
                                        Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1, voter);
                                    end

                                    break % both sides are already filled out!
                                else
                                    ihatvals=vals{1,side}.*fibsval{1,side}(usedidx{voter,side},patientsel);
                                    Ihat(test,side,voter) = ea_discfibers_getpeak(ihatvals(test), obj.posvisible, obj.negvisible, 'peak5');

                                    Ihat_train_global(numTestIt,training,side,voter) = ea_discfibers_getpeak(ihatvals(training), obj.posvisible, obj.negvisible, 'peak5');

                                    if isstruct(obj.ADJ)
                                        disp('Adjacency matrix for lateral symptoms is currently not supported')
                                    end
                                end
                        end
                    case {'Correlations / E-fields (Irmen 2020)', 'Reverse T-Tests / E-Fields (binary vars)', 'Odds Ratios / EF-Sigmoid (Jergas 2023)','Weighted Linear Regression / EF-Sigmoid (Dembek 2023)'} % efields
                        switch lower(obj.basepredictionon)
                            case 'profile of scores: spearman'
                                if lateral_score == false
                                    Ihat_all = corr(vals_flat,fibsval_usedidx_flat,'rows','pairwise','type','spearman');
                                    Ihat(test,1, voter) = Ihat_all(test);
                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_all(training);

                                    if isstruct(obj.ADJ)
                                        fibsval_ADJ_right = adj_scaler*obj.ADJ.ADJ(global_val_ind,global_conn_ind)*fibsval{voter,1}(1:end,patientsel);
                                        fibsval_ADJ_left = adj_scaler*obj.ADJ.ADJ(global_val_ind_lh,global_conn_ind_lh)*fibsval{voter,2}(1:end,patientsel);
                                        fibsval_ADJ_comb = [fibsval_ADJ_right;fibsval_ADJ_left];
                                        fibsval_ADJ_comb_full = fibsval_usedidx_flat + fibsval_ADJ_comb;

                                        Ihat(test,1,voter) = corr(vals_flat,fibsval_ADJ_comb_full(1:end,patientsel(test)),'rows','pairwise','type','spearman');
                                        Ihat(test,2, voter) = Ihat(test,1, voter);
                                        Ihat_train_global(numTestIt,training,1, voter) = corr(vals_flat,fibsval_ADJ_comb_full(1:end,patientsel(training)),'rows','pairwise','type','spearman');
                                        Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1, voter);
                                    end

                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = corr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(test)),'rows','pairwise','type','spearman');
                                    Ihat_train_global(numTestIt,training,side,voter) = corr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(training)),'rows','pairwise','type','spearman');
                                    if isstruct(obj.ADJ)
                                        disp('Adjacency matrix for lateral symptoms is currently not supported')
                                    end
                                end
                            case 'profile of scores: pearson'
                                if lateral_score == false
                                    Ihat_all = corr(vals_flat,fibsval_usedidx_flat,'rows','pairwise','type','pearson');
                                    Ihat(test,1, voter) = Ihat_all(test);
                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);

                                    if isstruct(obj.ADJ)
                                        fibsval_ADJ_right = adj_scaler*obj.ADJ.ADJ(global_val_ind,global_conn_ind)*fibsval{voter,1}(1:end,patientsel);
                                        fibsval_ADJ_left = adj_scaler*obj.ADJ.ADJ(global_val_ind_lh,global_conn_ind_lh)*fibsval{voter,2}(1:end,patientsel);
                                        fibsval_ADJ_comb = [fibsval_ADJ_right;fibsval_ADJ_left];
                                        fibsval_ADJ_comb_full = fibsval_usedidx_flat + fibsval_ADJ_comb;

                                        Ihat(test,1,voter) = corr(vals_flat,fibsval_ADJ_comb_full(1:end,patientsel(test)),'rows','pairwise','type','pearson');
                                        Ihat(test,2, voter) = Ihat(test,1, voter);
                                        Ihat_train_global(numTestIt,training,1, voter) = corr(vals_flat,fibsval_ADJ_comb_full(1:end,patientsel(training)),'rows','pairwise','type','pearson');
                                        Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1, voter);
                                    end

                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = corr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(test)),'rows','pairwise','type','pearson');
                                    Ihat_train_global(numTestIt,training,side,voter) = corr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(training)),'rows','pairwise','type','pearson');
                                    if isstruct(obj.ADJ)
                                        disp('Adjacency matrix for lateral symptoms is currently not supported')
                                    end
                                end
                            case 'profile of scores: bend'
                                if lateral_score == false
                                    Ihat_all = ea_bendcorr(vals_flat,fibsval_usedidx_flat);
                                    Ihat(test,1, voter) = Ihat_all(test);
                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);

                                    if isstruct(obj.ADJ)
                                        fibsval_ADJ_right = adj_scaler*obj.ADJ.ADJ(global_val_ind,global_conn_ind)*fibsval{voter,1}(1:end,patientsel);
                                        fibsval_ADJ_left = adj_scaler*obj.ADJ.ADJ(global_val_ind_lh,global_conn_ind_lh)*fibsval{voter,2}(1:end,patientsel);
                                        fibsval_ADJ_comb = [fibsval_ADJ_right;fibsval_ADJ_left];
                                        fibsval_ADJ_comb_full = fibsval_usedidx_flat + fibsval_ADJ_comb;

                                        Ihat(test,1,voter) = ea_bendcorr(vals_flat,fibsval_ADJ_comb_full(1:end,patientsel(test)));
                                        Ihat(test,2, voter) = Ihat(test,1, voter);
                                        Ihat_train_global(numTestIt,training,1, voter) = ea_bendcorr(vals_flat,fibsval_ADJ_comb_full(1:end,patientsel(training)));
                                        Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1, voter);
                                    end

                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_bendcorr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(test)));
                                    Ihat_train_global(numTestIt,training,side,voter) = ea_bendcorr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(training)));
                                    if isstruct(obj.ADJ)
                                        disp('Adjacency matrix for lateral symptoms is currently not supported')
                                    end
                                end
                            case 'mean of scores'
                                if lateral_score == false
                                    Ihat_all = ea_nanmean(vals_flat.*fibsval_usedidx_flat,1);
                                    Ihat(test,1, voter) = Ihat_all(test);
                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_all(training);
                                    if isstruct(obj.ADJ)
                                        Ihat_ADJ_rh = vals{voter,1}.*(obj.ADJ.ADJ(global_val_ind,global_conn_ind)*fibsval{1,1}(1:end,patientsel));
                                        Ihat_ADJ_lh = vals{voter,2}.*(obj.ADJ.ADJ(global_val_ind_lh,global_conn_ind_lh)*fibsval{1,2}(1:end,patientsel));
                                        Ihat_ADJ_comb = [Ihat_ADJ_rh;Ihat_ADJ_lh];
                                        Ihat_ADJ_test = transpose(adj_scaler*ea_nanmean(Ihat_ADJ_comb(1:end,test),1));
                                        Ihat_ADJ_training = adj_scaler*ea_nanmean(Ihat_ADJ_comb(1:end,training),1);
                                        Ihat(test,1,voter) = Ihat(test,1,voter) + Ihat_ADJ_test;
                                        Ihat(test,2, voter) = Ihat(test,1, voter);
                                        Ihat_train_global(numTestIt,training,1, voter) = Ihat_train_global(numTestIt,training,1, voter) + Ihat_ADJ_training;
                                        Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1, voter);
                                    end

                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                    Ihat_train_global(numTestIt,training,side,voter) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                    if isstruct(obj.ADJ)
                                        disp('Adjacency matrix for lateral symptoms is currently not supported')
                                    end
                                end
                            case 'sum of scores'
                                if lateral_score == false
                                    Ihat_all = ea_nansum(vals_flat.*fibsval_usedidx_flat,1);
                                    Ihat(test,1, voter) = Ihat_all(test);
                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);

                                    if isstruct(obj.ADJ)
                                        Ihat_ADJ_rh = vals{voter,1}.*(obj.ADJ.ADJ(global_val_ind,global_conn_ind)*fibsval{1,1}(1:end,patientsel));
                                        Ihat_ADJ_lh = vals{voter,2}.*(obj.ADJ.ADJ(global_val_ind_lh,global_conn_ind_lh)*fibsval{1,2}(1:end,patientsel));
                                        Ihat_ADJ_comb = [Ihat_ADJ_rh;Ihat_ADJ_lh];
                                        Ihat_ADJ_test = transpose(adj_scaler*ea_nansum(Ihat_ADJ_comb(1:end,test),1));
                                        Ihat_ADJ_training = adj_scaler*ea_nansum(Ihat_ADJ_comb(1:end,training),1);
                                        Ihat(test,1,voter) = Ihat(test,1,voter) + Ihat_ADJ_test;
                                        Ihat(test,2, voter) = Ihat(test,1, voter);
                                        Ihat_train_global(numTestIt,training,1, voter) = Ihat_train_global(numTestIt,training,1, voter) + Ihat_ADJ_training;
                                        Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1, voter);
                                    end

                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                    Ihat_train_global(numTestIt,training,side,voter) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                    if isstruct(obj.ADJ)
                                        disp('Adjacency matrix for lateral symptoms is currently not supported')
                                    end
                                end
                            case 'peak of scores'
                                if lateral_score == false
                                    Ihat_all = ea_discfibers_getpeak(vals_flat.*fibsval_usedidx_flat, obj.posvisible, obj.negvisible, 'peak');
                                    Ihat(test,1, voter) = Ihat_all(test);
                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);

                                    if isstruct(obj.ADJ)
                                        Ihat_ADJ_rh = vals{voter,1}.*(obj.ADJ.ADJ(global_val_ind,global_conn_ind)*fibsval{1,1}(1:end,patientsel));
                                        Ihat_ADJ_lh = vals{voter,2}.*(obj.ADJ.ADJ(global_val_ind_lh,global_conn_ind_lh)*fibsval{1,2}(1:end,patientsel));
                                        Ihat_ADJ_comb = [Ihat_ADJ_rh;Ihat_ADJ_lh];
                                        Ihat_ADJ_test = transpose(adj_scaler*ea_discfibers_getpeak(Ihat_ADJ_comb(1:end,test), obj.posvisible, obj.negvisible, 'peak'));
                                        Ihat_ADJ_training = adj_scaler*ea_discfibers_getpeak(Ihat_ADJ_comb(1:end,training), obj.posvisible, obj.negvisible, 'peak');
                                        Ihat(test,1,voter) = Ihat(test,1,voter) + Ihat_ADJ_test;
                                        Ihat(test,2, voter) = Ihat(test,1, voter);
                                        Ihat_train_global(numTestIt,training,1, voter) = Ihat_train_global(numTestIt,training,1, voter) + Ihat_ADJ_training;
                                        Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1, voter);
                                    end

                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_discfibers_getpeak(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)), obj.posvisible, obj.negvisible, 'peak');
                                    Ihat_train_global(numTestIt,training,side,voter) = ea_discfibers_getpeak(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)), obj.posvisible, obj.negvisible, 'peak');
                                    if isstruct(obj.ADJ)
                                        disp('Adjacency matrix for lateral symptoms is currently not supported')
                                    end
                                end
                            case 'peak 5% of scores'
                                if lateral_score == false
                                    Ihat_all = ea_discfibers_getpeak(vals_flat.*fibsval_usedidx_flat, obj.posvisible, obj.negvisible, 'peak5');
                                    Ihat(test,1, voter) = Ihat_all(test);
                                    Ihat(test,2, voter) = Ihat(test,1, voter);

                                    Ihat_train_global(numTestIt,training,1,voter) = Ihat_all(training);
                                    Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1,voter);

                                    if isstruct(obj.ADJ)
                                        Ihat_ADJ_rh = vals{voter,1}.*(obj.ADJ.ADJ(global_val_ind,global_conn_ind)*fibsval{1,1}(1:end,patientsel));
                                        Ihat_ADJ_lh = vals{voter,2}.*(obj.ADJ.ADJ(global_val_ind_lh,global_conn_ind_lh)*fibsval{1,2}(1:end,patientsel));
                                        Ihat_ADJ_comb = [Ihat_ADJ_rh;Ihat_ADJ_lh];

                                        Ihat_ADJ_test = transpose(adj_scaler*ea_discfibers_getpeak(Ihat_ADJ_comb(test), obj.posvisible, obj.negvisible, 'peak5'));
                                        Ihat_ADJ_training = adj_scaler*ea_discfibers_getpeak(Ihat_ADJ_comb(training), obj.posvisible, obj.negvisible, 'peak5');
                                        Ihat(test,1,voter) = Ihat(test,1,voter) + Ihat_ADJ_test;
                                        Ihat(test,2, voter) = Ihat(test,1, voter);
                                        Ihat_train_global(numTestIt,training,1, voter) = Ihat_train_global(numTestIt,training,1, voter) + Ihat_ADJ_training;
                                        Ihat_train_global(numTestIt,training,2, voter) = Ihat_train_global(numTestIt,training,1, voter);
                                    end

                                    break % both sides are already filled out!
                                else
                                    ihatvals=vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel);
                                    Ihat(test,side,voter) = ea_discfibers_getpeak(ihatvals(test), obj.posvisible, obj.negvisible, 'peak5');

                                    Ihat_train_global(numTestIt,training,side,voter) = ea_discfibers_getpeak(ihatvals(training), obj.posvisible, obj.negvisible, 'peak5');

                                    if isstruct(obj.ADJ)
                                        disp('Adjacency matrix for lateral symptoms is currently not supported')
                                    end
                                end
                        end
                end
            end
        end

    end
    %send out improvements of subscores
    for voter=1:size(vals,1)
        switch obj.multitractmode
            case 'Split & Color By Subscore'
                useI=obj.subscore.vars{voter};
            case 'Split & Color By PCA'
                useI=obj.subscore.pcavars{voter};
            otherwise
                useI=obj.responsevar;
        end
        for side=1:2
            mdl=fitglm(Ihat_train_global(numTestIt,training,side,voter),useI(training),lower(obj.predictionmodel));
            actualimprovs{voter,side}=predict(mdl,Ihat(test,side,voter));
        end
    end
    val_struct.vals=vals;
    val_struct.usedidx=usedidx;
try    val_struct.fibcell=fibcell; end % fibcell not always supplied.



end
