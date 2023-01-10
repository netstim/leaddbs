function [tractset]=ea_discfibers_optimize(tractset,app)
% Function to optimize parameters for fiber filtering using a Surrogate
% Optimization Algorithm. The app input is optional but can be used to
% live-view tuning of parameters.


%% The surrogate optimization algorithm alternates between two phases.

% Construct Surrogate — Create options.MinSurrogatePoints random points within the bounds. Evaluate the (expensive) objective function at these points. Construct a surrogate of the objective function by interpolating a radial basis function through these points.
% Search for Minimum — Search for a minimum of the objective function by sampling several thousand random points within the bounds. Evaluate a merit function based on the surrogate value at these points and on the distances between them and points where the (expensive) objective function has been evaluated. Choose the best point as a candidate, as measured by the merit function. Evaluate the objective function at the best candidate point. This point is called an adaptive point. Update the surrogate using this value and search again.


%% Definitions for Surrogate Optimization (from ML doc):
% The surrogate optimization algorithm description uses the following definitions.
%
% Current — The point where the objective function was evaluated most recently.
% Incumbent — The point with the smallest objective function value among all evaluated since the most recent surrogate reset.
% Best — The point with the smallest objective function value among all evaluated so far.
% Initial — The points, if any, that you pass to the solver in the InitialPoints option.
% Random points — Points in the Construct Surrogate phase where the solver evaluates the objective function. Generally, the solver takes these points from a quasirandom sequence, scaled and shifted to remain within the bounds. A quasirandom sequence is similar to a pseudorandom sequence such as rand returns, but is more evenly spaced. See https://en.wikipedia.org/wiki/Low-discrepancy_sequence. However, when the number of variables is above 500, the solver takes points from a Latin hypercube sequence. See https://en.wikipedia.org/wiki/Latin_hypercube_sampling.
% Adaptive points — Points in the Search for Minimum phase where the solver evaluates the objective function.
% Merit function — See Merit Function Definition.
% Evaluated points — All points at which the objective function value is known. These points include initial points, Construct Surrogate points, and Search for Minimum points at which the solver evaluates the objective function.
% Sample points. Pseudorandom points where the solver evaluates the merit function during the Search for Minimum phase. These points are not points at which the solver evaluates the objective function, except as described in Search for Minimum Details.

toolboxes_installed=ver;
if ~ismember('Global Optimization Toolbox',{toolboxes_installed.Name})
    ea_error('You need to install the Matlab Global Optimization Toolbox to use this functionality');
end

useparallel=0;

warning off
ks=[1,7,10]; %[2,7,10];
tractset.kIter=1; % make sure to run only 1 iteration.

% lower bounds for list of vars to optimize:
params=[1,  3,      1,  resolve_threshstrategy(tractset.threshstrategy)         % threshstrategy
    1,  3,      1,  resolve_corrtype(tractset.corrtype)                     % corrtype
    1,  4,      1,  resolve_efieldmetric(tractset.efieldmetric)             % efieldmetric
    1,  7,      1,  resolve_basepredictionon(tractset.basepredictionon)     % basepredictionon
    0,  100,    0,  mean(tractset.showposamount)                            % showposamount % will be multiplied by 10 for fixed amount
    0,  100,    0,  mean(tractset.shownegamount)                            % shownegamount % will be multiplied by 10 for fixed amount
    0,  100,    0,  tractset.connthreshold                                  % connthreshold
    1,  1,      0,  tractset.efieldthreshold];                              % efieldthreshold % will be scaled depending on efieldmetric

switch tractset.statmetric
    case {1,3,4,5} % t-tests, OSS-DBS, proportion tests, binomial tests
        paramidx=logical([1 % threshstrategy
            0 % corrtype
            0 % efieldmetric
            1 % basepredictionon
            1 % showposamount
            1 % shownegamount
            1 % connthreshold
            0]); % efieldthreshold
        params(4,2)=4; % upper bound change: baseprediction on does not allow profile of scores ideas for t-tests
        params(7,1:2)=[1,49]; % boundary change: connthreshold only allows 1-49% for t-tests
        params=params(paramidx,:);
    case {2,6,7} % e-fields, reverse-t-tests for binary & efields
        paramidx=logical([1 % threshstrategy
            1 % corrtype
            1 % efieldmetric
            1 % basepredictionon
            1 % showposamount
            1 % shownegamount
            1 % connthreshold
            1]); % efieldthreshold
        params=params(paramidx,:);
end
if tractset.statmetric==7
    choice=questdlg({'In many (but not all) cases, optimization in plain connections mode does not make too much sense, conceptually. Do you still want to proceed?'},'Plain Connections','Yes','No','No');
    if strcmp(choice,'No')
        return
    end
end


lb=params(:,1);
ub=params(:,2);
intcon=params(:,3);

% set up initial points with some good heuristics:
ip.X=params(:,4)'; % first initial point is whatever the user tried

ip.X=augmentips(ip.X);

%ip=repmat(ip',10,1);
%ip(2:end,1:4)=ip(2:end,1:4).*(1+(0.1*randn(size(ip(2:end,1:4)))));

options=optimoptions('surrogateopt',...
    'ObjectiveLimit',-0.9,... % optimal solution with average correlations of R~0.8 (rare to happen), lowest theoretical point is zero with an R of 1
    'MinSurrogatePoints',120,...
    'PlotFcn','surrogateoptplot',... 
    'Display','iter');
%    'CheckpointFile',fullfile(fileparts(tractset.leadgroup),'optimize_status.mat'),...

% check for parallel processing toolbox


if ismember('Parallel Computing Toolbox',{toolboxes_installed.Name}) && useparallel
    options.UseParallel=true;
    parpool('Processes',2);
end


% Solve problem
objconstr=@(x)struct('Fval',nestedfun(x));
if exist(fullfile(fileparts(tractset.leadgroup),'optimize_status.mat'),'file')
    choice=questdlg('Prior optimization has been done. Do you wish to continue on the same file?','Resume optimization?','Yes','Start from scratch','Yes');
    switch choice
        case 'Start from scratch'
        case 'Yes'
            priorstate=load(fullfile(fileparts(tractset.leadgroup),'optimize_status.mat'));
            ip=priorstate.ip;
    end
else
end

choice='y';
numIters=120;
while 1
    switch choice
        case 'y'
                options.InitialPoints=ip;
            if ~exist('numIters','var')
                numIters = input(sprintf('%s\n\n','Great, let us continue. How many trials do you want to run (enter amount)'),'s');
                numIters = str2double(numIters);
            end
            options.MaxFunctionEvaluations=numIters;
            [XOptim,fval,exitflag,output,ip]=surrogateopt(objconstr,lb,ub,find(intcon),options);
            save(fullfile(fileparts(tractset.leadgroup),'optimize_status.mat'),'ip');
        otherwise
            break
    end
    choice = input(sprintf('%s\n\n',['Optimal solution: Average R = ',num2str(-fval),'. Do you wish to continue optimizing? (y/n)']),'s');
    clear numIters
end
tractset=updatetractset(tractset,XOptim);



disp(['Optimal solution: Average R = ',num2str(-fval),'.']);
warning on

if ismember('Parallel Computing Toolbox',{toolboxes_installed.Name}) && useparallel
    poolobj = gcp('nocreate'); delete(poolobj);
end

tractset.save;

% Nested function that computes the objective function
    function Fval = nestedfun(X)

       tractset=updatetractset(tractset,X);

        if exist('app','var') % could be cool to update the state of the app GUI live.
            %updateapp(app,X);
        end

        Fval=getR(tractset);
    end

    function tractset=updatetractset(tractset,X)
        % resolve integer vars:
        disp('Parameters applied: ');

        tractset.threshstrategy=resolve_threshstrategy(X(1));
        tractset.corrtype=resolve_corrtype(X(2));
        tractset.efieldmetric=resolve_efieldmetric(X(3));
        tractset.basepredictionon=resolve_basepredictionon(X(4));

        fprintf('%s: %s\n','Threshold strategy',tractset.threshstrategy);
        fprintf('%s: %s\n','Correlation type',tractset.corrtype);
        fprintf('%s: %s\n','Efield metric',tractset.efieldmetric);
        fprintf('%s: %s\n','Base prediction on',tractset.basepredictionon);

        % set continuous vars:
        switch tractset.threshstrategy
            case 'Fixed Amount'
                tractset.showposamount=repmat(round(X(5))*10,1,2);
                tractset.shownegamount=repmat(round(X(6))*10,1,2);
            otherwise % percentages
                tractset.showposamount=repmat(X(5),1,2);
                tractset.shownegamount=repmat(X(6),1,2);
        end
        fprintf('%s: %01.1f \n','Show Positive Amount',tractset.showposamount);
        fprintf('%s: %01.1f \n','Show Negative Amount',tractset.shownegamount);


        tractset.connthreshold=X(7);
        fprintf('%s: %01.1f \n','Connectivity Threshold',tractset.connthreshold);

        switch tractset.efieldmetric % depending on efieldmetric, efieldthresholds should be scaled - input goes from 0 to 1
            case 'Sum'
                offset = 0.5;
                maxval = 5;
            case 'Mean'
                offset = 2.5;
                maxval = 1000;
            case {'Peak','Peak 5%'}
                offset = 0.05;
                maxval = 5;
        end
        tractset.efieldthreshold=(X(8)+offset)/(offset+1)*maxval;
        fprintf('%s: %01.1f \n','E-Field Threshold',tractset.efieldthreshold);

    end

    function R=getR(tractset)

        tractset.customselection = [];
        tractset.useExternalModel = false;
        for k=1:length(ks)
            if ks(k)==1 % circular
                tractset.customselection=tractset.patientselection;
                cvp.training{1}=ones(1,length(tractset.customselection));
                cvp.test{1}=ones(1,length(tractset.customselection));
                cvp.NumTestSets=1;
                try
                    [I, Ihat]=app.tractset.crossval(cvp,[],0,1);
                    if iscell(I)
                        for entry=1:length(I)
                            Rsub(entry)=corr(I{entry},Ihat{entry},'rows','pairwise');
                        end
                        R(R<0)=-1; % anything negative is the same in cross-validations
                        R(k)=ea_nanmean(Rsub);
                    else
                        R(k)=corr(I,Ihat,'rows','pairwise');
                    end
                catch
                    R(k)=nan;
                end
            else
                tractset.customselection=[];
                tractset.kfold = ks(k);
                try
                    [I,Ihat]=tractset.kfoldcv(1); % 1 = silent mode
                    if iscell(I)
                        for entry=1:length(I)
                            Rsub(entry)=corr(I{entry},Ihat{entry},'rows','pairwise');
                        end
                        R(R<0)=-1; % anything negative is the same in cross-validations
                        R(k)=ea_nanmean(Rsub);
                    else
                        R(k)=corr(I,Ihat,'rows','pairwise');
                    end
                catch
                    R(k)=nan;
                end
            end
        end
        R(R<0)=-1; % anything negative is the same in cross-validations
        R=ea_nanmean(R);
        if isnan(R) % out of bound settings
            R=-1; % minimal possible value
        end
        fprintf('%s: %01.0f\n','Average Correlation:',R);

        R=-R; % finally, flip, since we are minimizing. / could think instead to do R=1/exp(R) but less readible.

%         % map [-1:1] to a suitable gradient that downweights anything below 0 
%         % (since cross-validations are used):
%         zeropoint=1/(1+exp(5*(1-0.5)));
%         R=(1/(1+exp(5*(R-0.5))))-zeropoint;
    end





    function ip=augmentips(ip)
        switch tractset.statmetric
            case {1,3,4,5} % t-tests, OSS-DBS, proportion tests, binomial tests
                %% general heuristic 1: 20 % coverage
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    2 % basepredictionon: mean of scores
                    500 % showposamount: 500 on each side = 1000 fibers
                    500 % shownegamount: 500 on each side = 1000 fibers
                    20]' % connthreshold: 20% connected sites
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    2 % basepredictionon: mean of scores
                    1000 % showposamount: 1000 on each side = 2000 fibers
                    1000 % shownegamount: 1000 on each side = 2000 fibers
                    20]' % connthreshold: 20% connected sites
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    2 % basepredictionon: mean of scores
                    1000 % showposamount: 1000 on each side = 2000 fibers
                    200 % shownegamount: 200 on each side = 400 fibers
                    20]' % connthreshold: 20% connected sites
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    2 % basepredictionon: mean of scores
                    1000 % showposamount: 1000 on each side = 2000 fibers
                    0 % shownegamount: 0 on each side = 0 fibers
                    20]' % connthreshold: 20% connected sites
                    ];
                %% general heuristic 2: 30 % coverage
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    2 % basepredictionon: mean of scores
                    500 % showposamount: 500 on each side = 1000 fibers
                    500 % shownegamount: 500 on each side = 1000 fibers
                    30]' % connthreshold: 30% connected sites
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    2 % basepredictionon: mean of scores
                    1000 % showposamount: 1000 on each side = 2000 fibers
                    1000 % shownegamount: 1000 on each side = 2000 fibers
                    30]' % connthreshold: 30% connected sites
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    2 % basepredictionon: mean of scores
                    1000 % showposamount: 1000 on each side = 2000 fibers
                    200 % shownegamount: 200 on each side = 400 fibers
                    30]' % connthreshold: 30% connected sites
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    2 % basepredictionon: mean of scores
                    1000 % showposamount: 1000 on each side = 2000 fibers
                    0 % shownegamount: 0 on each side = 0 fibers
                    30]' % connthreshold: 30% connected sites
                    ];
            case {2,6,7} % e-fields, reverse-t-tests for binary & efields
                %% general heuristic 1: high tract impact, low % coverage
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    1 % corrtype: pearson
                    3 % efieldmetric: peak
                    3 % basepredictionon: peak of scores
                    500 % showposamount: 500 on each side = 1000 fibers
                    500 % shownegamount: 500 on each side = 1000 fibers
                    1 % connthreshold: 1% connected sites (at least one efield)
                    1]' % efieldthreshold: Peak: 1 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    1 % corrtype: pearson
                    3 % efieldmetric: peak
                    3 % basepredictionon: peak of scores
                    500 % showposamount: 500 on each side = 1000 fibers
                    250 % shownegamount: 250 on each side = 500 fibers
                    1 % connthreshold: 1% connected sites (at least one efield)
                    1]' % efieldthreshold: Peak: 1 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    1 % corrtype: pearson
                    3 % efieldmetric: peak
                    3 % basepredictionon: peak of scores
                    500 % showposamount: 500 on each side = 1000 fibers
                    0 % shownegamount: 0 on each side = 0 fibers
                    1 % connthreshold: 1% connected sites (at least one efield)
                    1]' % efieldthreshold: Peak: 1 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    1 % corrtype: pearson
                    3 % efieldmetric: peak
                    3 % basepredictionon: peak of scores
                    1000 % showposamount: 1000 on each side = 2000 fibers
                    0 % shownegamount: 0 on each side = 0 fibers
                    1 % connthreshold: 1% connected sites (at least one efield)
                    1]' % efieldthreshold: Peak: 1 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    1 % corrtype: pearson
                    3 % efieldmetric: peak
                    5 % basepredictionon: profile of scores: pearson
                    2000 % showposamount: 1000 on each side = 2000 fibers
                    2000 % shownegamount: 2000 on each side = 2000 fibers
                    1 % connthreshold: 1% connected sites (at least one efield)
                    1]' % efieldthreshold: Peak: 1 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    2 % corrtype: spearman
                    3 % efieldmetric: peak
                    7 % basepredictionon: profile of scores: spearman
                    2000 % showposamount: 1000 on each side = 2000 fibers
                    2000 % shownegamount: 2000 on each side = 2000 fibers
                    1 % connthreshold: 1% connected sites (at least one efield)
                    1]' % efieldthreshold: Peak: 1 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    3 % corrtype: bend
                    3 % efieldmetric: peak
                    7 % basepredictionon: profile of scores: bend
                    2000 % showposamount: 1000 on each side = 2000 fibers
                    2000 % shownegamount: 2000 on each side = 2000 fibers
                    1 % connthreshold: 1% connected sites (at least one efield)
                    1]' % efieldthreshold: Peak: 1 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    2 % corrtype: spearman
                    3 % efieldmetric: peak
                    3 % basepredictionon: peak of scores
                    500 % showposamount: 500 on each side = 1000 fibers
                    0 % shownegamount: 0 on each side = 0 fibers
                    1 % connthreshold: 1% connected sites (at least one efield)
                    1]' % efieldthreshold: Peak: 1 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    2 % corrtype: spearman
                    3 % efieldmetric: peak
                    3 % basepredictionon: peak of scores
                    500 % showposamount: 500 on each side = 1000 fibers
                    500 % shownegamount: 500 on each side = 1000 fibers
                    1 % connthreshold: 1% connected sites (at least one efield)
                    1]' % efieldthreshold: Peak: 1 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    3 % corrtype: bend
                    3 % efieldmetric: peak
                    3 % basepredictionon: peak of scores
                    500 % showposamount: 500 on each side = 1000 fibers
                    0 % shownegamount: 0 on each side = 0 fibers
                    1 % connthreshold: 1% connected sites (at least one efield)
                    1]' % efieldthreshold: Peak: 1 V/mm efieldthreshold
                    ];
                %% general heuristic 2: medium tract impact, 20 % coverage
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    1 % corrtype: pearson
                    3 % efieldmetric: peak
                    3 % basepredictionon: peak of scores
                    500 % showposamount: 500 on each side = 1000 fibers
                    500 % shownegamount: 500 on each side = 1000 fibers
                    20 % connthreshold: 20% connected sites
                    0.5]' % efieldthreshold: Peak: 0.5 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    1 % corrtype: pearson
                    3 % efieldmetric: peak
                    3 % basepredictionon: peak of scores
                    500 % showposamount: 500 on each side = 1000 fibers
                    250 % shownegamount: 250 on each side = 500 fibers
                    20 % connthreshold: 20% connected sites
                    0.5]' % efieldthreshold: Peak: 0.5 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    1 % corrtype: pearson
                    3 % efieldmetric: peak
                    3 % basepredictionon: peak of scores
                    500 % showposamount: 500 on each side = 1000 fibers
                    0 % shownegamount: 0 on each side = 0 fibers
                    20 % connthreshold: 20% connected sites
                    0.5]' % efieldthreshold: Peak: 0.5 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    1 % corrtype: pearson
                    3 % efieldmetric: peak
                    3 % basepredictionon: peak of scores
                    1000 % showposamount: 1000 on each side = 2000 fibers
                    0 % shownegamount: 0 on each side = 0 fibers
                    20 % connthreshold: 20% connected sites
                    0.5]' % efieldthreshold: Peak: 0.5 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    1 % corrtype: pearson
                    3 % efieldmetric: peak
                    5 % basepredictionon: profile of scores: pearson
                    2000 % showposamount: 1000 on each side = 2000 fibers
                    2000 % shownegamount: 2000 on each side = 2000 fibers
                    20 % connthreshold: 20% connected sites
                    0.5]' % efieldthreshold: Peak: 0.5 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    2 % corrtype: spearman
                    3 % efieldmetric: peak
                    7 % basepredictionon: profile of scores: spearman
                    2000 % showposamount: 1000 on each side = 2000 fibers
                    2000 % shownegamount: 2000 on each side = 2000 fibers
                    20 % connthreshold: 20% connected sites
                    0.5]' % efieldthreshold: Peak: 0.5 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    3 % corrtype: bend
                    3 % efieldmetric: peak
                    7 % basepredictionon: profile of scores: bend
                    2000 % showposamount: 1000 on each side = 2000 fibers
                    2000 % shownegamount: 2000 on each side = 2000 fibers
                    20 % connthreshold: 20% connected sites
                    0.5]' % efieldthreshold: Peak: 0.5 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    2 % corrtype: spearman
                    3 % efieldmetric: peak
                    3 % basepredictionon: peak of scores
                    500 % showposamount: 500 on each side = 1000 fibers
                    0 % shownegamount: 0 on each side = 0 fibers
                    20 % connthreshold: 20% connected sites
                    0.5]' % efieldthreshold: Peak: 0.5 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    2 % corrtype: spearman
                    3 % efieldmetric: peak
                    3 % basepredictionon: peak of scores
                    500 % showposamount: 500 on each side = 1000 fibers
                    500 % shownegamount: 500 on each side = 1000 fibers
                    20 % connthreshold: 20% connected sites
                    0.5]' % efieldthreshold: Peak: 0.5 V/mm efieldthreshold
                    ];
                %% variation of the above
                ip=[ip
                    [3 % thresh strategy: fixed amount
                    3 % corrtype: bend
                    3 % efieldmetric: peak
                    3 % basepredictionon: peak of scores
                    500 % showposamount: 500 on each side = 1000 fibers
                    0 % shownegamount: 0 on each side = 0 fibers
                    20 % connthreshold: 20% connected sites
                    0.5]' % efieldthreshold: Peak: 0.5 V/mm efieldthreshold
                    ];
        end

    end



    function out=resolve_efieldmetric(out)
        if ischar(out)
            switch out
                case 'Sum'
                    out=1;
                case 'Mean'
                    out=2;
                case 'Peak'
                    out=3;
                case 'Peak 5%'
                    out=4;
            end
        else
            switch out
                case 1
                    out='Sum';
                case 2
                    out='Mean';
                case 3
                    out='Peak';
                case 4
                    out='Peak 5%';
            end
        end

    end

    function out=resolve_basepredictionon(out)
        if ischar(out)
            switch lower(out)
                case 'mean of scores'
                    out=1;
                case 'sum of scores'
                    out=2;
                case 'peak of scores'
                    out=3;
                case 'peak 5% of scores'
                    out=4;
                case 'profile of scores: spearman'
                    out=5;
                case 'profile of scores: pearson'
                    out=6;
                case 'profile of scores: bend'
                    out=7;
            end
        else
            switch out
                case 1
                    out='mean of scores';
                case 2
                    out='sum of scores';
                case 3
                    out='peak of scores';
                case 4
                    out='peak 5% of scores';
                case 5
                    out='profile of scores: spearman';
                case 6
                    out='profile of scores: pearson';
                case 7
                    out='profile of scores: bend';
            end
        end

    end

    function out=resolve_corrtype(out)
        if ischar(out)
            switch out
                case 'Pearson'
                    out=1;
                case 'Spearman'
                    out=2;
                case 'Bend'
                    out=3;
            end
        else
            switch out
                case 1
                    out='Pearson';
                case 2
                    out='Spearman';
                case 3
                    out='Bend';
            end
        end

    end
    function out=resolve_threshstrategy(out)
        if ischar(out)
            switch out
                case 'Percentage Relative to Peak'
                    out=1;
                case 'Percentage Relative to Amount'
                    out=2;
                case 'Fixed Amount'
                    out=3;
                case 'Histogram (CDF)'
                    out=4;
                case 'Fixed Fiber Value'
                    out=5;
            end
        else
            switch out
                case 1
                    out='Percentage Relative to Peak';
                case 2
                    out='Percentage Relative to Amount';
                case 3
                    out='Fixed Amount';
                case 4
                    out='Histogram (CDF)';
                case 5
                    out='Fixed Fiber Value';
            end
        end
    end


end

