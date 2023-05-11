function [tractsetclone]=ea_cleartune_optimize_surrogate(tractsetclone,patlist,app,~,~)
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

[paramsR,paramsL] = genParams(app);

lbR = [paramsR(:,1)];
lbL = [paramsL(:,1)];
ubR = [paramsR(:,2)];
ubL = [paramsL(:,2)];
intergCond=[paramsR(:,3)];
%max 8 contacts activated

A = [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0;0 0 0 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 -1 -1 0];
b = [8;-1];
Aeq = [];
beq = [];

% set up initial points with some good heuristics:
startptsR = [paramsR(:,4)];
startptsL = [paramsL(:,4)];% first initial point is activation in k2, with voltage of 0.5
optionsR=optimoptions('surrogateopt',...
    'ObjectiveLimit',-0.9,... % optimal solution with average Ihat ~0.9, lowest theoretical point is zero with an R of 1
    'MinSurrogatePoints',500,...
    'PlotFcn','surrogateoptplot',...
    'InitialPoints',startptsR,...
    'MaxFunctionEvaluations',10,...
    'Display','iter');
optionsL=optimoptions('surrogateopt',...
    'ObjectiveLimit',-0.9,... % optimal solution with average Ihat ~0.9, lowest theoretical point is zero with an R of 1
    'MinSurrogatePoints',500,...
    'PlotFcn','surrogateoptplot',...
    'InitialPoints',startptsL,...
    'MaxFunctionEvaluations',10,...
    'Display','iter');

%    'CheckpointFile',fullfile(fileparts(tractset.leadgroup),'optimize_status.mat'),...


% check for parallel processing toolbox
if ismember('Parallel Computing Toolbox',{toolboxes_installed.Name}) && useparallel
    optionsR.UseParallel=true;
    optionsL.UseParallel=true;
    parpool('Processes',2);
end


% Solve problem
objconstrR=@(x)struct('Fval',nestedfunR(app,x,patlist,1));
objconstrL=@(x)struct('Fval',nestedfunL(app,x,patlist,2));
%
choice='y';
numIters=10;
while 1
    switch choice
        case 'y'
            %options.InitialPoints=ip;
            if ~exist('numIters','var')
                numIters = input(sprintf('%s\n\n','Great, let us continue. How many trials do you want to run (enter amount)'),'s');
                numIters = str2double(numIters);
            end
            [XOptimR,fvalR,~,~,ipR]=surrogateopt(objconstrR,lbR,ubR,find(intergCond),A,b,Aeq,beq,optionsR);
            [XOptimL,fvalL,~,~,ipL]=surrogateopt(objconstrL,lbL,ubL,find(intergCond),A,b,Aeq,beq,optionsL);
            writeVTA = 1;
            XOptimR_reform = reformatX(XOptimR);
            XOptimL_reform = reformatX(XOptimL);
            inputsR = {patlist{1},XOptimR_reform(1),XOptimR_reform(2:end),0,1,writeVTA};
            inputsL = {patlist{1},XOptimL_reform(1),XOptimL_reform(2:end),0,2,writeVTA};
            ea_generate_optim_vat(inputsR{:});
            ea_generate_optim_vat(inputsL{:});
            XOptim = [XOptimR_reform,XOptimL_reform];
            for i=1:length(ipL)
                ip.X(i,:) = [ipR.X(i,:),ipL.X(i,:)];
                ip.FvalR(i,1) = ipR.Fval(i,1);
                ip.FvalL(i,1) = ipL.Fval(i,1);
            end
            save(fullfile(patlist{1},'optimize_status_surrogate.mat'),'ip','XOptim');
        otherwise
            break
    end
    choice = input(sprintf('%s\n\n',['Optimal predicted Ihat R = ',num2str(-1*fvalR),' Optimal predicted Ihat L = ',num2str(-1*fvalL),'.Do you wish to continue optimizing? (y/n)']),'s');
    clear numIters
end
%tractsetclone=updatetractset(tractsetclone,XOptim);
avgIhat = ((-1*fvalR)+(-1*fvalL))/2;
disp(['Optimal solution: Average Ihat(R,L) = ',num2str(avgIhat),'.']);
warning on
if ismember('Parallel Computing Toolbox',{toolboxes_installed.Name}) && useparallel
    poolobj = gcp('nocreate'); delete(poolobj);
end

tractsetclone.save;

function [paramsR,paramsL] = genParams(app)

perc_contactdist = repmat([-100, 100, 0, 0],app.inputVars.numContacts,1);
bool_contact = repmat([0,1,1,0],app.inputVars.numContacts,1);
perc_contactdist(app.inputVars.startContact,end) = 100;
bool_contact(app.inputVars.startContact,end) = 1;
case_contact = [0,1,1,1];
amplitude = [app.inputVars.minCurr, app.inputVars.maxCurr, 0, app.inputVars.minCurr];
paramsR = [amplitude;perc_contactdist;bool_contact;case_contact];
paramsL = paramsR;

return
end

% Nested function that computes the objective function
function Fval = nestedfunR(app,X,patlist,side)
X = reformatX(X);
if isempty(X)
    Fval = Inf;
    return
end
constCurr = 0;
tractsetclone=updateStim(tractsetclone,X);
Fval=getFval(app,X,patlist,constCurr,tractsetclone,side);
return
end

function Fval = nestedfunL(app,X,patlist,side)
X = reformatX(X);
if isempty(X)
    Fval = Inf;
    return
end
constCurr = 0;
tractsetclone=updateStim(tractsetclone,X);
Fval=getFval(app,X,patlist,constCurr,tractsetclone,side);
return
end

function X = reformatX(X)
%if the indicator is off, the contact should also be off
for xx=10:17
    if X(xx) == 0
        X(xx-8) = 0;
    end
end
for xx=2:9
    if abs(X(xx)) < 10
        X(xx) = 0;
        X(xx+8) = 0;
    end
    
end

% Method 2 - scaling vars
if all(~X(2:9)) %all contacts have zero % activation its not allowed
    X = [];
    return
end

%doing it this way because harcoding indices is wrong and won't give us
%flexibility for future
%right
whichpos = X > 0;
whichpos(1) = 0;
whichpos(10:end) = 0;
whichneg = X < 0;
whichneg(1) = 0;
whichneg(10:end) = 0;
if any(whichpos)
    sum_whichpos_r = sum(X(whichpos));
    X(whichpos) = (X(whichpos)*100)/(sum_whichpos_r);
    X(end) = 0;
else
    X(end) = 100;
end
if any(whichneg)
    sum_whichneg_r = sum(X(whichneg));
    X(whichneg) = (X(whichneg)*-100)/(sum_whichneg_r);
end

return
end

function tractsetclone=updateStim(tractsetclone,X)
    disp('Parameters applied: ');
    fprintf('Amplitude:%d\n',X(1));
    whichContact = find(X(10:17));
    fprintf('Active contact: k0%d\n',whichContact-1)
    fprintf('Percentage activation: %d\n',X(2:9));
    fprintf('Case activation: %d\n',X(end));
end

function Fval=getFval(app,X,patlist,constCurr,tractsetclone,side)
    
    %create a vta inside this function, send it to cleartune, get Ihat out and return it as
    %Fval
    for yy=1:length(patlist)
        ampsel = X(1);
        concsel = [X(2:end)];
        writeVTA = 0;
        modelVTA = app.inputVars.modelVTA;
        inputs = {patlist{yy},ampsel,concsel,constCurr,side,writeVTA,modelVTA};
        [Efields,allS]=ea_generate_optim_vat(inputs{:});
        app.protocol{yy}.inputs=inputs;
        app.protocol{yy}.Efields=Efields;
        app.protocol{yy}.allS=allS;
        tractsetclone=ea_disctract;
        app.tractset.copyobj(tractsetclone);
    end
    [~,Ihat,actualimprovs] = runcrossval(app,'suggest',tractsetclone,patlist,side);
    Ihat=[Ihat(1:length(Ihat)/2),Ihat(length(Ihat)/2+1:end)]; % Ihat is exported as a column vector for both sides. Reformat to Nx2.
    if sum(isnan(Ihat(:)))/length(Ihat(:))>0.3
        ea_warning(['Many (',num2str(sum(isnan(Ihat(:)))*100/length(Ihat(:))),' percent) stimulation settings were not covered well by the model. Stimulation suggestions may not be meaningful. Please adjust model parameters in ',...
            app.fibfiltmodelpath,'.']);
    end
    preFval = calculateFval(tractsetclone,Ihat,actualimprovs,side);
    Fval = -1*preFval(side);

    
end

end

function preFval = calculateFval(tractsetclone,Ihat,actualimprovs,side)
    weightmatrix=zeros(size(actualimprovs,1),1); % in cleartune case always the same weights for any side and "patient" (which is VTA)
    for voter=1:length(weightmatrix)
         % same weight for all subjects in that voter (slider was used)
            weightmatrix(voter)=tractsetclone.cleartunevars.weights(1,voter);
    end
    weightmatrix_sum = ea_nansum(weightmatrix);
    for xx=1:size(weightmatrix,1) % make sure voter weights sum up to 1
        % for xx=1:size(Ihat_voters,1) % make sure voter weights sum up to 1
        % for yy=1:size(Ihat_voters,2)
        weightmatrix(xx)=weightmatrix(xx)./weightmatrix_sum;
    end
    for i=1:size(Ihat,1)
        for j=1:size(Ihat,2)
            wt_Ihat(i,j) = Ihat(i,j).*weightmatrix(i);
        end
    end
    preFval = ea_nansum(wt_Ihat,side); %should be the same since we are doing only one side now
    return
end

%     for i=1:size(actualimprovs,1)
%         for j=1:size(actualimprovs,2)
%             wt_actualimprovs{i,j} = actualimprovs{i,j}.*weightmatrix(i);
%         end
%     end
%     for side = 1:2
%         wtavg_actualimprovs(:,side) = ea_nansum([wt_actualimprovs{1:size(wt_actualimprovs,1),side}],2);
%     end
%     for side=1:2
%        maxImprv(side)=ea_nanmax(wtavg_actualimprovs(:,side));
%     end
%    preFval = maxImprv;
%    return

%end


