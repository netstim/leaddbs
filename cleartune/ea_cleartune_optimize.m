function [tractsetclone]=ea_cleartune_optimize(tractsetclone,patlist,app,opt_file,electrode_type)
% Function to optimize parameters for cleartune using a Surrogate
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

% if ~exist('command','var') || isempty(command)
%     command = 'cv';
% end
% 
% % When loading prior optimization, load and stop by default
% if ~exist('mode','var') || isempty(mode)
%     mode = 'stop';
% end

toolboxes_installed=ver;
if ~ismember('Global Optimization Toolbox',{toolboxes_installed.Name})
    ea_error('You need to install the Matlab Global Optimization Toolbox to use this functionality');
end

useparallel=0;

warning off
%app = ea_cleartune;
% list of vars to optimize:



params = [0.5,  5, 0, 1  % voltage, currently its only in mA
            -100, 100, 0, 0 %activation in k0
            -100, 100, 0, 0 %activation in k1
            -100, 100, 0, 0 %activation in k2
            -100, 100, 0, -100 %activation in k3
            -100, 100, 0,  0 %activation in k4
            -100, 100, 0,  0 %activation in k5
            -100, 100, 0,  0 %activation in k6
            -100, 100, 0,  0 %activation in k7
             0, 1, 1, 0 %k0 0 is off; 1 is pos; 2 is neg
             0, 1, 1, 0 %k1
             0, 1, 1, 0 %k2
             0, 1, 1, 1 %k3
             0, 1, 1, 0 %k4
             0, 1, 1, 0 %k5
             0, 1, 1, 0 %k6
             0, 1, 1, 0 %k7
             0, 1, 1, 1]; %case


% params = [0.5,  5, 0, 0.5  % voltage, currently its only in mA
%             0, 100, 0, 0 %activation in k0
%             0, 100, 0, 0 %activation in k1
%             0, 100, 0, 0 %activation in k2
%             -100, 0, 0, -100 %activation in k3
%             -100, 0, 0,  0 %activation in k4
%             -100, 0, 0,  0 %activation in k5
%             -100, 0, 0,  0 %activation in k6
%             -100, 0, 0,  0 %activation in k7
%             0, 1, 1, 1]; %case



lb=[params(:,1),params(:,1)];
ub=[params(:,2),params(:,2)];
intergCond=[params(:,3),params(:,3)];
%max 8 contacts activated

A = [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0;0 0 0 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 -1 -1 0];
b = [8;-1];
Aeq = [];
beq = [];

% Aeq = [0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0];
% beq = 0;
% set up initial points with some good heuristics:
startpts=[params(:,4),params(:,4)]; % first initial point is activation in k2, with voltage of 0.5
options=optimoptions('surrogateopt',...
    'ObjectiveLimit',-0.9,... % optimal solution with average Ihat ~0.9, lowest theoretical point is zero with an R of 1
    'MinSurrogatePoints',200,...
    'PlotFcn','surrogateoptplot',...
    'InitialPoints',startpts ,...
    'Display','iter');
%    'CheckpointFile',fullfile(fileparts(tractset.leadgroup),'optimize_status.mat'),...


% check for parallel processing toolbox
if ismember('Parallel Computing Toolbox',{toolboxes_installed.Name}) && useparallel
    options.UseParallel=true;
    parpool('Processes',2);
end

% Solve problem
objconstr=@(x)struct('Fval',nestedfun(x,patlist));
%
choice='y';
numIters=200;
while 1
    switch choice
        case 'y'
            %options.InitialPoints=ip;
            if ~exist('numIters','var')
                numIters = input(sprintf('%s\n\n','Great, let us continue. How many trials do you want to run (enter amount)'),'s');
                numIters = str2double(numIters);
            end
            options.MaxFunctionEvaluations=numIters;
            options.MinSampleDistance = 0.05;
            [XOptim,fval,exitflag,output,ip]=surrogateopt(objconstr,lb,ub,find(intergCond),A,b,Aeq,beq,options);
            
            save(fullfile(fileparts(tractsetclone.leadgroup),'optimize_status.mat'),'ip');
        otherwise
            break
    end
    choice = input(sprintf('%s\n\n',['Optimal predicted Ihat = ',num2str(fval) '.Do you wish to continue optimizing? (y/n)']),'s');
    clear numIters
end
%tractsetclone=updatetractset(tractsetclone,XOptim);

disp(['Optimal solution: Average Ihat(R,L) = ',num2str(fval),'.']);

warning on
if ismember('Parallel Computing Toolbox',{toolboxes_installed.Name}) && useparallel
    poolobj = gcp('nocreate'); delete(poolobj);
end

%tractsetclone.save;

% Nested function that computes the objective function
function Fval = nestedfun(X,patlist)
X_right = X(1:length(X)/2);
X_left = X(length(X)/2+1:end);

%if the indicator is off, the contact should also be off
for i=10:17
    if X_right(i) == 0
        X_right(i-8) = 0;
    end
    if X_left(i) == 0
        X_left(i-8) = 0;
    end
end

% Method 2 - scaling vars
if all(~X_right(2:9)) || all(~X_left(2:9))%all contacts have zero % activation its not allowed
    Fval = Inf;
    return
end
%doing it this way because harcoding indices is wrong and won't give us
%flexibility for future
%right
whichpos_r = X_right > 0;
whichpos_r(1) = 0;
whichpos_r(10:end) = 0;
whichneg_r = X_right < 0;
whichneg_r(1) = 0;
whichneg_r(10:end) = 0;
if any(whichpos_r)
    sum_whichpos_r = sum(X_right(whichpos_r));
    X_right(whichpos_r) = (X_right(whichpos_r)*100)/(sum_whichpos_r);
    X_right(end) = 0;
else
    X_right(end) = 100;
end
if any(whichneg_r)
    sum_whichneg_r = sum(X_right(whichneg_r));
    X_right(whichneg_r) = (X_right(whichneg_r)*-100)/(sum_whichneg_r);
end
%left
whichpos_l = X_left > 0;
whichpos_l(1) = 0;
whichpos_l(10:end) = 0;
whichneg_l = X_left < 0;
whichneg_l(1) = 0;
whichneg_l(10:end) = 0;
if any(whichpos_l)
    sum_whichpos_l = sum(X_left(whichpos_l));
    X_left(whichpos_l) = (X_left(whichpos_l)*100)/(sum_whichpos_l);
    X_left(end) = 0;
else
    X_left(end) = 100;
end
if any(whichneg_l)
    sum_whichneg_l = sum(X_left(whichneg_l));
    X_left(whichneg_l) = (X_left(whichneg_l)*-100)/(sum_whichneg_l);
end


constCurr = 0;
tractsetclone=updateStim(tractsetclone,X_left,X_right);
Fval=getFval(X_left,X_right,patlist,constCurr,tractsetclone);
return
end

function tractsetclone=updateStim(tractsetclone,X_left,X_right)
    disp('Parameters applied R: ');
    fprintf('Amplitude R:%d\n',X_right(1));
    fprintf('Amplitude L:%d\n',X_left(1));
    whichContactR = find(X_right(10:17));
    whichContactL = find(X_left(10:17));
    fprintf('Active contact R: k0%d\n',whichContactR-1)
    fprintf('Active contact L: k0%d\n',whichContactL-1)
    fprintf('Percentage activation R: %d\n',X_right(2:10));
    fprintf('Percentage activation L: %d\n',X_left(2:10));

end

function Fval=getFval(X_left,X_right,patlist,constCurr,tractsetclone)
    
    %create a vta inside this function, send it to cleartune, get Ihat out and return it as
    %Fval
    
    for i=1:length(patlist)
        ampsel_r = X_right(1);
        ampsel_l = X_left(1);
        concsel_r = [X_right(2:9),X_right(end)];
        concsel_l = [X_left(2:9),X_left(end)];
        inputs = {patlist{i},ampsel_r,ampsel_l,concsel_r,concsel_l,constCurr};
        [Efields,allS]=ea_generate_optim_vat(inputs{:});
        app.protocol{i}.inputs=inputs;
        app.protocol{i}.Efields=Efields;
        app.protocol{i}.allS=allS;
        tractsetclone=ea_disctract;
        app.tractset.copyobj(tractsetclone);
    end
    [I,Ihat,actualimprovs] = runcrossval(app,'suggest',tractsetclone,patlist);
    preFval = calculateFval(tractsetclone,Ihat,actualimprovs);
    disp(preFval);
    Fval = -1*mean(preFval);

    
end

end

function preFval = calculateFval(tractsetclone,Ihat,actualimprovs)
    weightmatrix=zeros(size(actualimprovs,1),1); % in cleartune case always the same weights for any side and "patient" (which is VTA)
    for voter=1:length(weightmatrix)
         % same weight for all subjects in that voter (slider was used)
            weightmatrix(voter)=tractsetclone.cleartunevars.weights(voter);
    end
    weightmatrix_sum = ea_nansum(weightmatrix);
    for xx=1:size(weightmatrix,1) % make sure voter weights sum up to 1
        % for xx=1:size(Ihat_voters,1) % make sure voter weights sum up to 1
        % for yy=1:size(Ihat_voters,2)
        weightmatrix(xx)=weightmatrix(xx)./weightmatrix_sum;
    end
    for i=1:size(actualimprovs,1)
        for j=1:size(actualimprovs,2)
            wt_actualimprovs{i,j} = actualimprovs{i,j}.*weightmatrix(i);
        end
    end
    for side = 1:2
        wtavg_actualimprovs(:,side) = ea_nansum([wt_actualimprovs{1:size(wt_actualimprovs,1),side}],2);
        
    end
    for side=1:2
       maxImprv(side)=ea_nanmax(wtavg_actualimprovs(:,side));
    end
    preFval = maxImprv;
    return
end


