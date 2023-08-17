function ea_save_fibscore_model(obj, ExternalModelFile)

    % this function is called from the GUI

    % we want to define the model outside of cross-validation
    if ~exist('patsel','var') % patsel can be supplied directly (in this case, obj.patientselection is ignored), e.g. for cross-validations.
        patientsel = obj.patientselection;
    end

    % get fiber model (vals) and corresponding indices (usedidx)
    if obj.cvlivevisualize
        [vals,fibcell,usedidx] = ea_discfibers_calcstats(obj, patientsel);
        obj.draw(vals,fibcell,usedidx)
        drawnow;
    else
        [vals,~,usedidx] = ea_discfibers_calcstats(obj, patientsel);
    end

    % maps fiber model (vals) into the full connectome space
    vals_all = cell(size(vals)); 
    for voter=1:size(vals,1)
        for side=1:size(vals,2)
            vals_all{voter,side} = zeros(obj.results.(ea_conn2connid(obj.connectome)).totalFibers,1);
            switch obj.connectivity_type
                case 2
                    vals_all{voter,side}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{side}(usedidx{voter,side})) = vals{voter,side};
                otherwise
                    vals_all{voter,side}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{side}(usedidx{voter,side})) = vals{voter,side};
            end
        end
    end

    % get Ihats and fit the linear model
    % here we compute Ihats for all patietsel at once
    disp('Computing the linear model...')
    fibsval = full(obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval);

    switch obj.multitractmode
        case {'Single Tract Analysis', 'Split & Color By Subscore'}
            % Ihat is the estimate of improvements (not scaled to real improvements)
            if strcmp(obj.multitractmode,'Single Tract Analysis')
                Ihat = nan(length(patientsel),2);
                Ihat_ph = nan(length(patientsel),2);
            else
                Ihat = nan(length(patientsel),2,length(obj.subscore.vars));
                Ihat_ph = nan(length(patientsel),2,length(obj.subscore.vars));
            end
            training = 1:size(patientsel,2);
            test = training;
            Improvement = obj.responsevar(patientsel,:);
            [Ihat, ~, ~,~] = ea_compute_fibscore_model(1, obj.adj_scaler, obj, fibsval, Ihat, Ihat_ph, patientsel, training, test);
            predictor=squeeze(ea_nanmean(Ihat,2));
            mdl=fitglm(predictor(training),Improvement(training),lower(obj.predictionmodel));
            Intercept = mdl.Coefficients.Estimate(1);
            Slope = mdl.Coefficients.Estimate(2);
            plotName = 'Fitting of a linear model for the stored Fiber Model';
            empiricallabel = 'Variable of Intrest';
            pred_label = 'Fiber (Ihat) score';    
            LM_values_slope = ['Slope: ' num2str(Slope)];
            LM_values_intercept = ['Intercept: ' num2str(Intercept)];
            ftr.mdl = mdl;
        
            h=ea_corrbox(Improvement,predictor,'permutation',{plotName,empiricallabel,pred_label, plotName, LM_values_slope, LM_values_intercept});
    
        case 'Split & Color By PCA'
            Ihat    = nan(length(patientsel),2,length(obj.subscore.vars));
            Ihat_ph = nan(length(patientsel),2,length(obj.subscore.vars));
            training = 1:size(patientsel,2);
            test = training;
            for i=1:length(obj.subscore.vars)
                Improvement{i} = obj.subscore.vars{i}(patientsel);
            end

            [Ihat, ~, ~,~] = ea_compute_fibscore_model(1, obj.adj_scaler, obj, fibsval, Ihat, Ihat_ph, patientsel, training, test);
            % NB : here we have just predictor because training and test are the same
            predictor=squeeze(ea_nanmean(Ihat,2));

            for voter=1:size(obj.subscore.pcavars,1)
                 useI=obj.subscore.pcavars{voter}(patientsel);
                 covariates=[];
                 for cv = 1:length(obj.covars)
                     covariates = [covariates,obj.covars{cv}(patientsel)];
                 end
                 if ~isempty(covariates)
                     mdl=fitglm([predictor(training,voter),covariates(training,:)],useI(training),lower(obj.predictionmodel));
                 else
                     mdl=fitglm([predictor(training,voter)],useI(training),lower(obj.predictionmodel));
                 end
                 % save model 
                 ftr.mdl{voter} = mdl;
                 if ~isempty(covariates)
                     Ihat_ph(test,:,voter)=repmat(predict(mdl,[predictor(test,voter),covariates(test,:)]),1,2); % fill both sides equally
                 else
                     Ihat_ph(test,:,voter)=repmat(predict(mdl,[predictor(test,voter)]),1,2); % fill both sides equally
                 end
            end
            
            % Quick check : plot PC scores prediction
            Ihat=Ihat_ph; 
            Ihat=squeeze(ea_nanmean(Ihat,2));
            for pcc=1:obj.numpcs
                if obj.subscore.posvisible(pcc)==1 || obj.subscore.negvisible(pcc)==1 % don't try to plot if not showing any fibers for this PC
                    ea_corrplot(obj.subscore.pcavars{pcc}(patientsel),Ihat(:,pcc), 'noperm', ...
                        {['Disc. Fiber prediction for PC ',num2str(pcc)],'PC score (Empirical)','PC score (Predicted)'},...
                        [], [], obj.subscore.pcacolors(pcc, :));
                end
            end

            % Save PCA weights, variable names, covariate names 
            ftr.pcacoeff = obj.subscore.pcacoeff; 
            ftr.pcavarnames = obj.subscore.labels; 
            ftr.covarnames = obj.covarlabels; 
            ftr.numpc = size(obj.subscore.pcavars,1); 

    end 

    % add the connectome name for recognition
    ftr.vals_all = vals_all;
    ftr.connectome = ea_conn2connid(obj.connectome);
    ftr.conn_type = obj.connectivity_type;
    ftr.fibsvalType = ea_method2methodid(obj);

    save(ExternalModelFile, '-struct', 'ftr');

end
