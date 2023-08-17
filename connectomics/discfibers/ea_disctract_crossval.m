function [I,Ihat,cvs,val_struct]=ea_disctract_crossval(viz,tractset,strategy,iterations,posthoccorrect,customconfig)

if ~exist('customconfig','var')
    customconfig=[];
end

switch tractset.multitractmode
    case 'Split & Color By Group'
        if ismember(strategy(1:6),{'Custom','Import'})
            ea_error('Multi Group cross-val not implemented for Custom selections or Imported Models.')
        end
        usedgroups=unique(tractset.M.patient.group(tractset.patientselection));
        allptsel=tractset.patientselection;
        tractset.multitractmode='Single Tract Analysis';
        try % this section is in try catch since dangerous if errors out (would then change patsel and multitractmode without the user knowing).
            for g=1:length(usedgroups)
                thisgrouppatients=find(tractset.M.patient.group==usedgroups(g));
                is=ismember(allptsel,thisgrouppatients);
                tractset.patientselection=allptsel(is);
                [I{g},Ihat{g},cvs,sel,val_struct{g}]=ea_disctract_crossval_do(~viz,tractset,strategy,iterations,customconfig);
                if viz
                    ea_disctract_crossval_visualize(tractset,I{g},Ihat{g},cvs,posthoccorrect,sel,usedgroups(g));
                end
            end
        catch
            disp('Multitract crossvalidation failed.');
        end
        tractset.multitractmode='Split & Color By Group';
        tractset.patientselection=allptsel;
    otherwise
        [I,Ihat,cvs,sel,val_struct{1}]=ea_disctract_crossval_do(~viz,tractset,strategy,iterations,customconfig);
        if viz
            ea_disctract_crossval_visualize(tractset,I,Ihat,cvs,posthoccorrect,sel);
        end
end



function [I,Ihat,cvs,sel,val_struct]=ea_disctract_crossval_do(silent,tractset,strategy,iterations,customconfig)
val_struct=nan;
switch strategy
    case 'Leave-Nothing-Out (Permutation-Based)'
        tractset.nestedLOO = false;
        tractset.useExternalModel = false;
        tractset.customselection = [];

        if ~isfield(customconfig, 'permcorrtype')
            customconfig.permcorrtype = 'Spearman';
        end
        [I,Ihat,R0,R1,pperm,~,val_struct]=tractset.lnopb(customconfig.permcorrtype,silent);
        if ~silent
            if strcmp(tractset.multitractmode,'Split & Color By PCA')
                I=mat2cell( squeeze(I(1,:,:)), length(tractset.patientselection), ones(1, length(tractset.subscore.vars)));
                Ihat=Ihat{1};
                for subvar = 1:length(tractset.subscore.vars)
                    h1=ea_plothistperm([tractset.subscore.labels{subvar},' [Permutation-Based Test]'],[R1(subvar);R0(:,subvar)], ...
                        {'Unpermuted prediction'},{1},1);
                    try saveas(h1,[fileparts(tractset.leadgroup),filesep,'fiberfiltering',filesep,tractset.ID,'_PCA_',tractset.subscore.labels{subvar},'_',cvs,'_permtest.png']); end
                end

            else
                I = I(:,1);
                Ihat = Ihat{1};
                h1=ea_plothistperm([tractset.responsevarlabel,' [Permutation-Based Test]'],[R1;R0],{'Unpermuted prediction'},{1},1);
                try saveas(h1,[fileparts(tractset.leadgroup),filesep,'fiberfiltering',filesep,tractset.ID,'_',tractset.responsevarlabel,'_',cvs,'_permtest.png']); end
            end
        end
        sel=tractset.patientselection;

        cvs='lnopb';


    case 'Leave-One-Patient-Out'
        tractset.customselection = [];
        tractset.useExternalModel = false;
        [I,Ihat]=tractset.loocv;
        cvs='loocv';
        sel=tractset.patientselection;
    case 'Leave-One-Cohort-Out'
        tractset.nestedLOO = false;
        tractset.useExternalModel = false;
        tractset.customselection = [];
        [I,Ihat]=tractset.lococv;
        cvs='lococv';
        sel=tractset.patientselection;
    case 'k-fold (randomized)'
        if tractset.kfold==1 % circular
            tractset.nestedLOO = false;
            tractset.useExternalModel = false;
            tractset.customselection=tractset.patientselection;
            cvp.training{1}=true(1,length(tractset.customselection));
            cvp.test{1}=true(1,length(tractset.customselection));
            cvp.NumTestSets=1;
            [I, Ihat]=tractset.crossval(cvp,[],0,silent);
            cvs='circular';
            sel=tractset.patientselection;
        else
            tractset.customselection = [];
            tractset.useExternalModel = false;
            tractset.kIter = iterations;
            [I,Ihat,val_struct]=tractset.kfoldcv(silent); % val_struct used for optimizer
            cvs='kfoldcv';
            sel=tractset.patientselection;
        end
    case 'Custom (Patients)'
        if strcmp(tractset.multitractmode, 'Split & Color By PCA')
            ea_error('Please use a different CV strategy for PCA')
        else
            tractset.nestedLOO = false;
            tractset.useExternalModel = false;
            cvp.NumTestSets = 1;
            % training and test indices from the items list
            training = find(ismember(customconfig.trainonitems,customconfig.trainonvalues));
            test = find(ismember(customconfig.predictonitems,customconfig.predictonvalues));
            % Patient selected based on the training and test indices
            tractset.customselection = unique([training, test]);
            % Construct cvp struct
            cvp.training{1} = ismember(tractset.customselection, training);
            cvp.test{1} = ismember(tractset.customselection, test);
            [I, Ihat]=tractset.crossval(cvp,[],0,silent);
            cvs = 'custom_pts';
            sel = test;
        end
    case 'Import Model'
        tractset.nestedLOO = false;
        cvp.NumTestSets = 1;
        % only test indices from the items list
        test = find(ismember(customconfig.predictonitems,customconfig.predictonvalues));
        % also use them for training as a place holder
        training = find(ismember(customconfig.predictonitems,customconfig.predictonvalues));
        tractset.useExternalModel = true;
        % Patient selected based on the training and test indices
        tractset.customselection = unique([training, test]);
        % Construct cvp struct
        cvp.training{1} = ismember(tractset.customselection, training);
        cvp.test{1} = ismember(tractset.customselection, test);
        [I, Ihat]=tractset.crossval(cvp);
        cvs = 'imported_model';
        sel = test;
    case 'Custom (Cohorts)'
        if strcmp(tractset.multitractmode, 'Split & Color By PCA')
            ea_error('Please use a different CV strategy for PCA')
        else
            tractset.nestedLOO = false;
            tractset.useExternalModel = false;
            cvp.NumTestSets = 1;
            % training and test indices from the items list
            training = find(ismember(tractset.M.patient.group, str2double(customconfig.trainonvalues)))';
            test = find(ismember(tractset.M.patient.group,str2double(customconfig.predictonvalues)))';
            % Patient selected based on the training and test indices
            tractset.customselection = unique([training, test]);
            % Construct cvp struct
            cvp.training{1} = ismember(tractset.customselection, training);
            cvp.test{1} = ismember(tractset.customselection, test);
            [I, Ihat]=tractset.crossval(cvp);
            cvs = 'custom_cohs';
            sel = test;
        end
    case 'Custom (Subcohorts)'
        tractset.nestedLOO = false;
        tractset.useExternalModel = false;
        cvp.NumTestSets = 1;
        % training and test indices from the items list
        tsets = find(ismember(tractset.setlabels, customconfig.trainonvalues));
        training=[];
        for ts=1:length(tsets)
            training = [training,find(tractset.setselections{tsets(ts)})];
        end
        training=unique(training);
        tsets = find(ismember(tractset.setlabels, customconfig.predictonvalues));
        test=[];
        for ts=1:length(tsets)
            test = [test,find(tractset.setselections{tsets(ts)})];
        end
        test=unique(test);
        % Patient selected based on the training and test indices
        tractset.customselection = unique([training, test]);
        % Construct cvp struct
        cvp.training{1} = ismember(tractset.customselection, training);
        cvp.test{1} = ismember(tractset.customselection, test);

        % When in PCA mode, training set should remain consistent with patient selection
        % because the PCA is computed on the patient selection
        if strcmp(tractset.multitractmode, 'Split & Color By PCA')
            if ~isequal(training, tractset.patientselection)
                ea_error(['You should train the model on the patients used to compute the PCA, ' ...
                    'otherwise the PCA coefficients will not correspond!'])
            end
        end

        [I, Ihat]=tractset.crossval(cvp);
        cvs = 'custom_subcohs';
        sel = test;
    case 'Custom (Sets)'
        if strcmp(tractset.multitractmode, 'Split & Color By PCA')
            ea_error('Please use a different CV strategy for PCA')
        else
            tractset.nestedLOO = false;
            tractset.useExternalModel = false;
            cvp.NumTestSets = 1;
            % Generate cvpartition
            rng(tractset.rngseed);
            cvpKfold = cvpartition(length(tractset.M.patient.list), 'KFold', str2double(kfoldk));
            % Calculate training and test indices based on the sets selected and the cvpartitions generated
            trainSetInd = str2double(customconfig.trainonvalues');
            testSetInd = str2double(customconfig.predictonvalues');
            training = zeros(size(tractset.M.patient.list));
            test = zeros(size(tractset.M.patient.list));
            for i=1:length(trainSetInd)
                training = training | cvpKfold.test(trainSetInd(i));
            end
            for i=1:length(testSetInd)
                test = test | cvpKfold.test(testSetInd(i));
            end
            training = find(training);
            test = find(test);
            % Patient selected based on the training and test indices
            tractset.customselection = unique([training, test]);
            % Construct cvp struct
            cvp.training{1} = ismember(tractset.customselection, training);
            cvp.test{1} = ismember(tractset.customselection, test);
            [I, Ihat]=tractset.crossval(cvp);

            cvs = 'custom_sets';
            sel = test;
        end
end