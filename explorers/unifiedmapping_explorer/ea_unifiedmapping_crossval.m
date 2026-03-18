function [I,Ihat,cvs,val_struct]=ea_unifiedmapping_crossval(viz,explorer,strategy,iterations,posthoccorrect,customconfig)

if ~exist('customconfig','var')
    customconfig=[];
end

switch explorer.multitractmode
    case 'Split & Color By Group'
        if ismember(strategy(1:6),{'Custom','Import'})
            ea_error('Multi Group cross-val not implemented for Custom selections or Imported Models.')
        end
        usedgroups=unique(explorer.M.patient.group(explorer.patientselection));
        allptsel=explorer.patientselection;
        explorer.multitractmode='Single Tract Analysis';
        try % this section is in try catch since dangerous if errors out (would then change patsel and multitractmode without the user knowing).
            for g=1:length(usedgroups)
                thisgrouppatients=find(explorer.M.patient.group==usedgroups(g));
                is=ismember(allptsel,thisgrouppatients);
                explorer.patientselection=allptsel(is);
                [I{g},Ihat{g},cvs,sel,val_struct{g}]=ea_unifiedmapping_crossval_do(~viz,explorer,strategy,iterations,customconfig);
                if viz
                    ea_unifiedmapping_crossval_visualize(explorer,I{g},Ihat{g},cvs,posthoccorrect,sel,usedgroups(g), toolUsed);
                end
            end
        catch
            disp('Multitract crossvalidation failed.');
        end
        explorer.multitractmode='Split & Color By Group';
        explorer.patientselection=allptsel;
    otherwise
        cnt=0;
        if strcmp(explorer.activated.sweetspotmapping,'On')
            cnt=cnt+1;
            explorer.drawTool='sweetspotmapping';
            [I{cnt},Ihat{cnt},cvs{cnt},sel{cnt},val_struct{cnt}{1}]=ea_unifiedmapping_crossval_do(~viz,explorer,strategy,iterations,customconfig);
            if viz
                toolUsed=1;
                ea_unifiedmapping_crossval_visualize(explorer,I{cnt},Ihat{cnt},cvs{cnt},posthoccorrect,sel{cnt}, [], toolUsed);
            end
        end
        if strcmp(explorer.activated.fiberfiltering,'On')
            cnt=cnt+1;
            explorer.drawTool='fiberfiltering';
            [I{cnt},Ihat{cnt},cvs{cnt},sel{cnt},val_struct{cnt}{1}]=ea_unifiedmapping_crossval_do(~viz,explorer,strategy,iterations,customconfig);
            if viz
                toolUsed = 2;
                ea_unifiedmapping_crossval_visualize(explorer,I{cnt},Ihat{cnt},cvs{cnt},posthoccorrect,sel{cnt}, [], toolUsed);
            end
        end
        if strcmp(explorer.activated.networkmapping,'On')
            cnt=cnt+1;
            explorer.drawTool='networkmapping';
            [I{cnt},Ihat{cnt},cvs{cnt},sel{cnt},val_struct{cnt}{1}]=ea_unifiedmapping_crossval_do(~viz,explorer,strategy,iterations,customconfig);
            if viz
                toolUsed = 3;
                ea_unifiedmapping_crossval_visualize(explorer,I{cnt},Ihat{cnt},cvs{cnt},posthoccorrect,sel{cnt}, [], toolUsed);
            end            
        end

        if cnt>1
            ea_glmplot(cell2mat(Ihat),I{1},'Normal',{'Unified Model','Estimates','Empirical'});
        end

end



function [I,Ihat,cvs,sel,val_struct]=ea_unifiedmapping_crossval_do(silent,explorer,strategy,iterations,customconfig)
val_struct=nan;
switch strategy
    case 'Leave-Nothing-Out (Permutation-Based)'
        explorer.nestedLOO = false;
        explorer.useExternalModel = false;
        explorer.customselection = [];

        if ~isfield(customconfig, 'permcorrtype')
            customconfig.permcorrtype = 'Spearman';
        end
        [I,Ihat,R0,R1,pperm,~,val_struct]=explorer.lnopb(customconfig.permcorrtype,silent);
        if ~silent
            if strcmp(explorer.multitractmode,'Split & Color By PCA')
                I=mat2cell( squeeze(I(1,:,:)), length(explorer.patientselection), ones(1, length(explorer.subscore.vars)));
                Ihat=Ihat{1};
                for subvar = 1:length(explorer.subscore.vars)
                    h1=ea_plothistperm([explorer.subscore.labels{subvar},' [Permutation-Based Test]'],[R1(subvar);R0(:,subvar)], ...
                        {'Unpermuted prediction'},{1},1);
                    try saveas(h1,[fileparts(explorer.leadgroup),filesep,'fiberfiltering',filesep,explorer.ID,'_PCA_',explorer.subscore.labels{subvar},'_',cvs,'_permtest.png']); end
                end

            else
                I = I(:,1);
                Ihat = Ihat{1};
                h1=ea_plothistperm([explorer.responsevarlabel,' [Permutation-Based Test]'],[R1;R0],{'Unpermuted prediction'},{1},1);
                try saveas(h1,[fileparts(explorer.leadgroup),filesep,'fiberfiltering',filesep,explorer.ID,'_',explorer.responsevarlabel,'_',cvs,'_permtest.png']); end
            end
        end
        sel=explorer.patientselection;

        cvs='lnopb';


    case 'Leave-One-Patient-Out'
        explorer.customselection = [];
        explorer.useExternalModel = false;
        [I,Ihat]=explorer.loocv;
        cvs='loocv';
        sel=explorer.patientselection;
    case 'Leave-One-Cohort-Out'
        explorer.nestedLOO = false;
        explorer.useExternalModel = false;
        explorer.customselection = [];
        [I,Ihat]=explorer.lococv;
        cvs='lococv';
        sel=explorer.patientselection;
    case 'k-fold (randomized)'
        if explorer.kfold==1 % circular
            explorer.nestedLOO = false;
            explorer.useExternalModel = false;
            explorer.customselection=explorer.patientselection;
            cvp.training{1}=true(1,length(explorer.customselection));
            cvp.test{1}=true(1,length(explorer.customselection));
            cvp.NumTestSets=1;
            [I, Ihat]=explorer.crossval(cvp,[],0,silent);
            cvs='circular';
            sel=explorer.patientselection;
        else
            explorer.customselection = [];
            explorer.useExternalModel = false;
            explorer.kIter = iterations;
            [I,Ihat,val_struct]=explorer.kfoldcv(silent); % val_struct used for optimizer
            cvs='kfoldcv';
            sel=explorer.patientselection;
        end
    case 'Custom (Patients)'
        if strcmp(explorer.multitractmode, 'Split & Color By PCA')
            ea_error('Please use a different CV strategy for PCA')
        else
            explorer.nestedLOO = false;
            explorer.useExternalModel = false;
            cvp.NumTestSets = 1;
            % training and test indices from the items list
            training = find(ismember(customconfig.trainonitems,customconfig.trainonvalues));
            test = find(ismember(customconfig.predictonitems,customconfig.predictonvalues));
            % Patient selected based on the training and test indices
            explorer.customselection = unique([training, test]);
            % Construct cvp struct
            cvp.training{1} = ismember(explorer.customselection, training);
            cvp.test{1} = ismember(explorer.customselection, test);
            [I, Ihat]=explorer.crossval(cvp,[],0,silent);
            cvs = 'custom_pts';
            sel = test;
        end
    case 'Import Model'
        explorer.nestedLOO = false;
        cvp.NumTestSets = 1;
        % only test indices from the items list
        test = find(ismember(customconfig.predictonitems,customconfig.predictonvalues));
        % also use them for training as a place holder
        training = find(ismember(customconfig.predictonitems,customconfig.predictonvalues));
        explorer.useExternalModel = true;
        % Patient selected based on the training and test indices
        explorer.customselection = unique([training, test]);
        % Construct cvp struct
        cvp.training{1} = ismember(explorer.customselection, training);
        cvp.test{1} = ismember(explorer.customselection, test);
        [I, Ihat]=explorer.crossval(cvp);
        cvs = 'imported_model';
        sel = test;
    case 'Custom (Cohorts)'
        if strcmp(explorer.multitractmode, 'Split & Color By PCA')
            ea_error('Please use a different CV strategy for PCA')
        else
            explorer.nestedLOO = false;
            explorer.useExternalModel = false;
            cvp.NumTestSets = 1;
            % training and test indices from the items list
            training = find(ismember(explorer.M.patient.group, str2double(customconfig.trainonvalues)))';
            test = find(ismember(explorer.M.patient.group,str2double(customconfig.predictonvalues)))';
            % Patient selected based on the training and test indices
            explorer.customselection = unique([training, test]);
            % Construct cvp struct
            cvp.training{1} = ismember(explorer.customselection, training);
            cvp.test{1} = ismember(explorer.customselection, test);
            [I, Ihat]=explorer.crossval(cvp);
            cvs = 'custom_cohs';
            sel = test;
        end
    case 'Custom (Subcohorts)'
        explorer.nestedLOO = false;
        explorer.useExternalModel = false;
        cvp.NumTestSets = 1;
        % training and test indices from the items list
        tsets = find(ismember(explorer.setlabels, customconfig.trainonvalues));
        training=[];
        for ts=1:length(tsets)
            training = [training,find(explorer.setselections{tsets(ts)})];
        end
        training=unique(training);
        tsets = find(ismember(explorer.setlabels, customconfig.predictonvalues));
        test=[];
        for ts=1:length(tsets)
            test = [test,find(explorer.setselections{tsets(ts)})];
        end
        test=unique(test);
        % Patient selected based on the training and test indices
        explorer.customselection = unique([training, test]);
        % Construct cvp struct
        cvp.training{1} = ismember(explorer.customselection, training);
        cvp.test{1} = ismember(explorer.customselection, test);

        % When in PCA mode, training set should remain consistent with patient selection
        % because the PCA is computed on the patient selection
        if strcmp(explorer.multitractmode, 'Split & Color By PCA')
            if ~isequal(training, explorer.patientselection)
                ea_error(['You should train the model on the patients used to compute the PCA, ' ...
                    'otherwise the PCA coefficients will not correspond!']);
            end
        end

        [I, Ihat]=explorer.crossval(cvp);
        cvs = 'custom_subcohs';
        sel = test;
    case 'Custom (Sets)'
        if strcmp(explorer.multitractmode, 'Split & Color By PCA')
            ea_error('Please use a different CV strategy for PCA')
        else
            explorer.nestedLOO = false;
            explorer.useExternalModel = false;
            cvp.NumTestSets = 1;
            % Generate cvpartition
            rng(explorer.rngseed);
            cvpKfold = cvpartition(length(explorer.M.patient.list), 'KFold', str2double(kfoldk));
            % Calculate training and test indices based on the sets selected and the cvpartitions generated
            trainSetInd = str2double(customconfig.trainonvalues');
            testSetInd = str2double(customconfig.predictonvalues');
            training = zeros(size(explorer.M.patient.list));
            test = zeros(size(explorer.M.patient.list));
            for i=1:length(trainSetInd)
                training = training | cvpKfold.test(trainSetInd(i));
            end
            for i=1:length(testSetInd)
                test = test | cvpKfold.test(testSetInd(i));
            end
            training = find(training);
            test = find(test);
            % Patient selected based on the training and test indices
            explorer.customselection = unique([training, test]);
            % Construct cvp struct
            cvp.training{1} = ismember(explorer.customselection, training);
            cvp.test{1} = ismember(explorer.customselection, test);
            [I, Ihat]=explorer.crossval(cvp);

            cvs = 'custom_sets';
            sel = test;
        end
end