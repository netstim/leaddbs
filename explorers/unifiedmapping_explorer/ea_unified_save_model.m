function ea_unified_save_model(obj, ExternalModelFile, ss_val, ff_val, nm_val)

    exportedAnything = false;

    % we want to define the model outside of cross-validation
    if ~exist('obj.explorer.patientselection','var') % patsel can be supplied directly (in this case, obj.patientselection is ignored), e.g. for cross-validations.
        patientsel = obj.explorer.patientselection;
    %         obj.customselection=obj.patientselection;
    end
    
    if ss_val
        exportedAnything = true;
        sides = fieldnames(obj.explorer.spotdrawn);
        for i = 1:length(sides)
            spot{i} = obj.explorer.spotdrawn.(sides{i});
            exportmodel.sweetspot.model_vals{i} = spot{i}.img(:);
            exportmodel.sweetspot.space = spot{i};

        end
             % get Ihats and fit the linear model
        % here we compute Ihats for all patietsel at once
        % disp('Computing sweetspot linear model...')
        % % Ihat is the estimate of improvements (not scaled to real improvements)
        % training = 1:size(patientsel,2);
        % test = training;
        % cvp.training = {training}; %training & test are the same (lno)
        % cvp.test = {test};
        % cvp.NumTestSets = 1;
        % Improvement = obj.explorer.responsevar(patientsel,:);
        % [~, Ihat] = crossval(obj.explorer, cvp);
        % predictor=squeeze(ea_nanmean(Ihat,2));
        % mdl=fitglm(predictor(training),Improvement(training),'linear');
        % exportmodel.sweetspot.mdl = mdl;

    end
    
    if ff_val
        exportedAnything = true;
        fields=fieldnames(obj.explorer.fiberdrawn);
        patientsel = obj.explorer.patientselection;
        for f = 1:length(fields)
            exportmodel.fiberfiltering.(fields{f}) = obj.explorer.fiberdrawn.(fields{f});
        end
        exportmodel.fiberfiltering.connectome = ea_conn2connid(obj.explorer.calcsettings.fibfilt_connectome);
        % model fitting not added yet
    end

    if nm_val
        exportedAnything = true;
        network = obj.explorer.networkdrawn;
        exportmodel.networkmapping.model_vals=network.img;
        exportmodel.networkmapping.connectome = ea_conn2connid(obj.explorer.calcsettings.netmap_connectome);
        %     % get Ihats and fit the linear model
        % % here we compute Ihats for all patietsel at once
        % disp('Computing networkmapping linear model...')
        % % Ihat is the estimate of improvements (not scaled to real improvements)
        % training = 1:size(patientsel,2);
        % test = training;
        % cvp.training = {training}; %training & test are the same (lno)
        % cvp.test = {test};
        % cvp.NumTestSets = 1;
        % Improvement = obj.explorer.responsevar(patientsel,:);
        % [~, Ihat] = crossval(obj.explorer, cvp);
        % predictor=squeeze(ea_nanmean(Ihat,2));
        % mdl=fitglm(predictor(training),Improvement(training),'linear');
        % exportmodel.networkmapping.mdl = mdl;
    
    
        % add the connectome name for recognition
    %     nm.vals_all = vals_all;
        % exportmodel.networkmapping.model_vals=vals;

    end

    if ~exportedAnything
        warning('ea_unified_save_model: NoModels', 'No models were selected for export. Nothing was saved.')
        return;
    end

     %export results of space
    exportmodel.stattest = obj.explorer.statsettings;
    exportmodel.corrtype = obj.explorer.corrtype;
    exportmodel.e_field_metric = obj.explorer.e_field_metric;
    exportmodel.mirrored_sides = obj.explorer.mirrorsides;
    exportmodel.patientselection=patientsel;
    exportmodel.calcsettings = obj.explorer.calcsettings;
    exportmodel.analysispath = obj.explorer.analysispath;

    save(ExternalModelFile, '-struct', 'exportmodel');

end