function [I,Ihat,cvs]=ea_disctract_crossval(tractset,strategy,iterations,kfoldk,trainonitems,trainonvalues,predictonitems,predictonvalues)

switch strategy
    case 'Leave-Nothing-Out (Permutation-Based)'
        tractset.nestedLOO = false;
        tractset.useExternalModel = false;
        tractset.customselection = [];
        [I,Ihat,R0,R1,pperm]=tractset.lnopb;

        if strcmp(tractset.multitractmode,'Split & Color By PCA')
            I=mat2cell( squeeze(I(1,:,:)), length(tractset.patientselection), ones(1, length(tractset.subscore.vars)));
            Ihat=Ihat{1};
            for subvar = 1:length(tractset.subscore.vars)
                cvs='lnopb';
                sel=tractset.patientselection;
                h1=ea_plothistperm([tractset.subscore.labels{subvar},' [Permutation-Based Test]'],[R1(subvar);R0(:,subvar)], ...
                    {'Unpermuted prediction'},{1},1);
                try saveas(h1,[fileparts(tractset.leadgroup),filesep,'fiberfiltering',filesep,tractset.ID,'_PCA_',tractset.subscore.labels{subvar},'_',cvs,'_permtest.png']); end
            end

        else
            I = I(:,1);
            Ihat = Ihat{1};
            cvs='lnopb';
            sel=tractset.patientselection;
            h1=ea_plothistperm([tractset.responsevarlabel,' [Permutation-Based Test]'],[R1;R0],{'Unpermuted prediction'},{1},1);
            try saveas(h1,[fileparts(tractset.leadgroup),filesep,'fiberfiltering',filesep,tractset.ID,'_',tractset.responsevarlabel,'_',cvs,'_permtest.png']); end
        end

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
        tractset.customselection = [];
        tractset.useExternalModel = false;
        tractset.kIter = iterations;
        [I,Ihat]=tractset.kfoldcv;
        cvs='kfoldcv';
        sel=tractset.patientselection;
    case 'Custom (Patients)'
        if strcmp(tractset.multitractmode, 'Split & Color By PCA')
            ea_error('Please use a different CV strategy for PCA')
        else
            tractset.nestedLOO = false;
            tractset.useExternalModel = false;
            cvp.NumTestSets = 1;
            % training and test indices from the items list
            training = find(ismember(trainonitems,trainonvalues));
            test = find(ismember(predictonitems,predictonvalues));
            % Patient selected based on the training and test indices
            tractset.customselection = unique([training, test]);
            % Construct cvp struct
            cvp.training{1} = ismember(tractset.customselection, training);
            cvp.test{1} = ismember(tractset.customselection, test);
            [I, Ihat]=tractset.crossval(cvp);
            cvs = 'custom_pts';
            sel = test;
        end
    case 'Import Model'
        tractset.nestedLOO = false;
        cvp.NumTestSets = 1;
        % only test indices from the items list
        test = find(ismember(predictonitems,predictonvalues));
        % also use them for training as a place holder
        training = find(ismember(predictonitems,predictonvalues));
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
            training = find(ismember(tractset.M.patient.group, str2double(trainonvalues)))';
            test = find(ismember(tractset.M.patient.group,str2double(predictonvalues)))';
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
        tsets = find(ismember(tractset.setlabels, trainonvalues));
        training=[];
        for ts=1:length(tsets)
            training = [training,find(tractset.setselections{tsets(ts)})];
        end
        training=unique(training);
        tsets = find(ismember(tractset.setlabels, predictonvalues));
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
            trainSetInd = str2double(trainonvalues');
            testSetInd = str2double(predictonvalues');
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