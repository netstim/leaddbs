function ea_save_networkmapping_model(obj, ExternalModelFile)

    % this function is called from the GUI

    % we want to define the model outside of cross-validation
    if ~exist('patsel','var') % patsel can be supplied directly (in this case, obj.patientselection is ignored), e.g. for cross-validations.
        patientsel = obj.patientselection;
    end

    % get the selected netmap model (vals) and corresponding indices (usedidx)
    if obj.cvlivevisualize
        vals = ea_networkmapping_calcstats(obj, patientsel);
        obj.draw(vals)
        drawnow;
    else
        [vals] = ea_networkmapping_calcstats(obj, patientsel);
    end

%     if obj.mirrorsides
%         getpatsel=[patientsel,patientsel+length(obj.allpatients)];
%     end

%     % get the raw fingerprints from the model
%     vals_all = obj.results.(ea_conn2connid(obj.connectome)).connval(getpatsel,:); 

    % get Ihats and fit the linear model
    % here we compute Ihats for all patietsel at once
    disp('Computing the linear model...')
    % Ihat is the estimate of improvements (not scaled to real improvements)
    training = 1:size(patientsel,2);
    test = training;
    cvp.training = {training}; %training & test are the same (lno)
    cvp.test = {test};
    cvp.NumTestSets = 1;
    Improvement = obj.responsevar(patientsel,:);
    [~, Ihat] = crossval(obj, cvp);
    predictor=squeeze(ea_nanmean(Ihat,2));
    mdl=fitglm(predictor(training),Improvement(training),'linear');
    Intercept = mdl.Coefficients.Estimate(1);
    Slope = mdl.Coefficients.Estimate(2);
    plotName = 'Fitting of a linear model for the stored Network Model';
    empiricallabel = 'Variable of Intrest';
    pred_label = 'Network (Ihat) score';
    LM_values_slope = ['Slope: ' num2str(Slope)];
    LM_values_intercept = ['Intercept: ' num2str(Intercept)];
    nm.mdl = mdl;


    % add the connectome name for recognition
%     nm.vals_all = vals_all;
    nm.model_vals=vals;
    nm.connectome = ea_conn2connid(obj.connectome);
    nm.model_type = obj.statmetric;
    nm.mirrored_sides = obj.mirrorsides;
    nm.corrtype = obj.corrtype;
    save(ExternalModelFile, '-struct', 'nm');


end