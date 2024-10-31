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
%list of vars to optimize:
% [paramsR,paramsL] = genParams(app);
%
% lbR = [paramsR(:,1)];
% lbL = [paramsL(:,1)];
% ubR = [paramsR(:,2)];
% ubL = [paramsL(:,2)];


%for bipolar
%max 8 contacts activated
%A = [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0;0 0 0 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 -1 -1 0];
%b = [8;-1];

% A = [];
% b = [];
% 8 mA if cathodic stimulation
% set up initial points with some good heuristics:
% startptsR = [paramsR(:,4)];
% startptsL = [paramsL(:,4)];% first initial point is activation in k2, with voltage of 0.5


% New way to defide the stimulation window
% we assume one patient
% check the reconstruction file (we are assuming the same model on both sides)

for pt = 1:length(patlist)
    ptindx = pt;
    %we need the reconstruction file for getting current per contact.
    %Defines the segmented and ring contacts
    reconstruction_file = dir([patlist{pt},filesep,'reconstruction',filesep,'sub-*desc-reconstruction.mat']);
    reconstruction_file_path = [reconstruction_file.folder,filesep,reconstruction_file.name];
    [reconst, ~, ~, ~] = ea_get_reconstruction(reconstruction_file_path);
    % do not update S here, just get the bounds in mA!
    [min_bound_per_contactR, max_bound_per_contactR, ~] = ea_get_currents_per_contact(app.MinCylindricEditField_2.Value,app.MaxCylindricEditField.Value, app.MinSegmentedEditField.Value, app.MaxSegmentedEditField.Value, reconst, 0, 0);
    min_bound_per_contactL = min_bound_per_contactR;
    max_bound_per_contactL = max_bound_per_contactR;
    %arbritary start points
    startcontactR = 5;
    startcontactL = 3;
    startptsR = zeros(1,app.inputVars.numContacts);
    startptsL = zeros(1,app.inputVars.numContacts);
    % set third contact (k2) to the middle of the higher current bound
    if abs(max_bound_per_contactR(startcontactR)) > abs(min_bound_per_contactR(startcontactR))
        startptsR(startcontactR) = max_bound_per_contactR(startcontactR)/1.0;
    else
        startptsR(startcontactR) = min_bound_per_contactR(startcontactR)/1.0;
    end

    if abs(max_bound_per_contactL(startcontactL)) > abs(min_bound_per_contactL(startcontactL))
        startptsL(startcontactL) = max_bound_per_contactL(startcontactL)/1.0;
    else
        startptsL(startcontactL) = min_bound_per_contactL(startcontactL)/1.0;
    end
    %define lower bounds and upper bounds
    lbR = min_bound_per_contactR;
    lbL = min_bound_per_contactL;
    ubR = max_bound_per_contactR;
    ubL = max_bound_per_contactL;
    %the current values per contact are not integers. They can assume
    %decimal values
    intergCond = zeros(1,app.inputVars.numContacts);
    %the number of contacts is automatically obtained from your
    %reconstruction file
    A = repmat(-1,1,app.inputVars.numContacts);
    b = 5;
    Aeq = [];
    beq = [];
    %define the right side & left side options
    optionsR=optimoptions('surrogateopt',...
        'ObjectiveLimit',-0.9,... % optimal solution with average Ihat ~0.9, lowest theoretical point is zero with an R of 1
        'MinSurrogatePoints',200,...
        'PlotFcn','surrogateoptplot',...
        'InitialPoints',startptsR,...
        'MaxFunctionEvaluations',500,... %choose 500 iters for now
        'Display','iter');
    optionsL=optimoptions('surrogateopt',...
        'ObjectiveLimit',-0.9,... % optimal solution with average Ihat ~0.9, lowest theoretical point is zero with an R of 1
        'MinSurrogatePoints',200,...
        'PlotFcn','surrogateoptplot',...
        'InitialPoints',startptsL,...
        'MaxFunctionEvaluations',500,... %choose 500 iters for now
        'Display','iter');

    % check for parallel processing toolbox
    if ismember('Parallel Computing Toolbox',{toolboxes_installed.Name}) && useparallel
        optionsR.UseParallel=true;
        optionsL.UseParallel=true;
        parpool('Processes',2);
    end

    % for monopolar case
    if all(ubR(:) <= 0.0)
        objconstrR=@(x)struct('Fval',nestedfunR_Monopolar(x,patlist,ptindx));
        [XOptimR,fvalR,~,~,ipR]=surrogateopt(objconstrR,lbR,ubR,find(intergCond),A,b,Aeq,beq,optionsR);
        if length(ipR.Fval) == 2
            keyboard
        end
    else %bipolar, not well tested
        objconstrR=@(x)struct('Fval',nestedfunR(x,patlist,ptindx));
        [XOptimR,fvalR,~,~,ipR]=surrogateopt(objconstrR,lbR,ubR,optionsR);
    end
    %to output VTA
    writeVTA = 1;
    %default: SimBio/FieldTrip
    modelVTA = app.inputVars.modelVTA;
    %remove very small currents
    XOptimR(abs(XOptimR) < 0.000001) = 0.0;
    %because the inputs are normalized to remove small currents
    newoptimR = reformatX(XOptimR);
    if ~isempty(newoptimR)
        XOptimR = newoptimR;
    end
    for i=1:size(ipR.X,1)
        newR = reformatX(ipR.X(i,:));
        if ~isempty(newR)
            ipR.X(i,:) = newR;
        end
    end
    save(fullfile(patlist{pt},'optimize_status_surrogate_R.mat'),'ipR','XOptimR');
    % generate the VTA
    [inputsR,ampl_R,perc_val_R] = ea_get_inputs_for_optimizer(patlist,XOptimR, modelVTA,writeVTA,1);
    ea_generate_optim_vat(inputsR{:});

    if all(ubL(:) <= 0.0)
        objconstrL=@(x)struct('Fval',nestedfunL_Monopolar(x,patlist,ptindx));
        [XOptimL,fvalL,~,~,ipL]=surrogateopt(objconstrL,lbL,ubL,find(intergCond),A,b,Aeq,beq,optionsL);
    else
        objconstrL=@(x)struct('Fval',nestedfunL(x,patlist,ptindx));
        [XOptimL,fvalL,~,~,ipL]=surrogateopt(objconstrL,lbL,ubL,optionsL);
    end

    XOptimL(abs(XOptimL) < 0.000001) = 0.0;

    newoptimL = reformatX(XOptimL);
    if ~isempty(newoptimL)
        XOptimL = newoptimL;
    end
    for i=1:size(ipL.X,1)
        newL = reformatX(ipL.X(i,:));
        if ~isempty(newL)
            ipL.X(i,:) = newL;
        end
    end
    save(fullfile(patlist{pt},'optimize_status_surrogate_L.mat'),'ipL','XOptimL');
    [inputsL,ampl_L,perc_val_L] = ea_get_inputs_for_optimizer(patlist{pt},XOptimL, app.inputVars.modelVTA,writeVTA,2);
    ea_generate_optim_vat(inputsL{:});



    %finally, we can save the stimulation parameters of the winning VTA
    options = setOPTS(patlist{pt});
    S = ea_initializeS(options);
    S = ea_cleartune_generateMfile([ampl_R,perc_val_R],[ampl_L,perc_val_L],S,0);
    save(fullfile(patlist{pt},'desc-stimparameters.mat'),'S');

    %create a plot that shows how Ihat changes with amplitude
    createIhatAmpPlot(app,inputsR,inputsL)
    avgIhat = ((-1*fvalR)+(-1*fvalL))/2;
    disp(['Optimal solution: Average Ihat(R,L) = ',num2str(avgIhat),'.']);
    warning on
    %end the parallel toolbox
    if ismember('Parallel Computing Toolbox',{toolboxes_installed.Name}) && useparallel
        poolobj = gcp('nocreate'); delete(poolobj);
    end

end

% Nested function that computes the objective function
function f = nestedfunR(X,patlist,ptindx)

% limit to 8mA
if sum(X(X>0)) > sum(abs(X(X<0)))
    f.Ineq = sum(X(X>0)) - 8.0;
else
    f.Ineq = sum(abs(X(X<0))) - 8.0;
end

X = reformatX(X);

if isempty(X)
    f.Fval = Inf;
    return
end
constCurr = 2;
tractsetclone=updateStim(tractsetclone,X);
f.Fval=getFval(app,X,patlist,1,ptindx);

return
end

function f = nestedfunL(X,patlist,ptindx)
X = reformatX(X);

if sum(X(X>0)) > sum(abs(X(X<0)))
    f.Ineq = sum(X(X>0)) - 8.0;
else
    f.Ineq = sum(abs(X(X<0))) - 8.0;
end

if isempty(X)
    f.Fval = Inf;
    return
end
constCurr = 2;
tractsetclone=updateStim(tractsetclone,X);
f.Fval=getFval(app,X,patlist,2,ptindx);
% limit to 8mA

return
end

% Nested function that computes the objective function
function Fval = nestedfunR_Monopolar(X,patlist,ptindx)

X = reformatX(X);

if isempty(X)
    Fval = Inf;
    return
end
tractsetclone=updateStim(tractsetclone,X);
Fval=getFval(app,X,patlist,1,ptindx);

return
end

function Fval = nestedfunL_Monopolar(X,patlist,ptindx)
X = reformatX(X);

if isempty(X)
    Fval = Inf;
    return
end

tractsetclone=updateStim(tractsetclone,X);
Fval=getFval(app,X,patlist,2,ptindx);
% limit to 8mA

return
end

function X = reformatX(X)
    if all(~X(1:8)) %all contacts have zero % activation its not allowed
        X = [];
        return
    end
return
end

function tractsetclone=updateStim(tractsetclone,X)
    disp('Parameters applied: ');
    %fprintf('Amplitude:%d\n',X(1));
    whichContact = find(X(1:end));
    fprintf('Active contact: k0%d\n',whichContact-1)
    fprintf('Current, mA: %d\n',X(1:end));
    %disp('Case activation:100%');
end

function Fval=getFval(app,X,patlist,side,ptindx)

    %create a vta inside this function, send it to cleartune, get Ihat out and return it as
    %Fval
    selpat = patlist{ptindx};
    % remove very small currents
    X(abs(X) < 0.000001) = 0.0;
    writeVTA = 0;

    % new way to define inputs
    inputs = ea_get_inputs_for_optimizer(selpat,X, app.inputVars.modelVTA,writeVTA,side);
    try
        [Efields,allS]=ea_generate_optim_vat(inputs{:});
    catch ME
        if (contains(ME.message,'Despite all attempts the VTA model could not be created'))
            %try once more
            [Efields,allS]=ea_generate_optim_vat(inputs{:});
            %Fval =  NaN;

        end
    end
    app.protocol{ptindx}.inputs=inputs;
    app.protocol{ptindx}.Efields=Efields;
    app.protocol{ptindx}.allS=allS;
    tractsetclone=ea_disctract;
    app.tractset.copyobj(tractsetclone);
    [~,Ihat,actualimprovs] = runcrossval(app,'suggest',tractsetclone,patlist,ptindx,side);
    %Ihat=[Ihat(1:length(Ihat)/2),Ihat(length(Ihat)/2+1:end)]; % Ihat is exported as a column vector for both sides. Reformat to Nx2.
    if isnan(Ihat(:))
        ea_warning("Some test values have returned NaNs..why?");
    end
    preFval = calculateFval(app,Ihat,actualimprovs,side,ptindx);
    Fval = -1*preFval;
    %add penalty function if user chooses
    Fval = penaltyFunc(app,X,Fval);
return

end

end

function preFval = calculateFval(app,Ihat,actualimprovs,side,ptindx)
    weightmatrix=zeros(size(actualimprovs,1),1); % in cleartune case always the same weights for any side and "patient" (which is VTA)
    for voter=1:length(weightmatrix)
         % same weight for all subjects in that voter (slider was used)
            weightmatrix(voter)=app.symptomWeightVar{ptindx,side}(voter);
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
    preFval = ea_nansum(wt_Ihat(:,side)); %should be the same since we are doing only one side now

    return
end
function Fval = penaltyFunc(app,X,Fval)
    if strcmp(app.ApplypenaltyusingDropDown.Value,'Quadratic curve')
        xx = app.SweetspotamplitudeEditField.Value;
        yy = abs(sum(X));
        Fval = Fval + ((yy-xx)^2)/100;
    elseif strcmp(app.ApplypenaltyusingDropDown.Value,'Range of amplitudes')
        if sum(abs(X)) > 3 || sum(abs(X)) < 2
            Fval = Fval + app.ApplyapenaltyvalueofEditField.Value;
        end
    else
        return
    end
end
function createIhatAmpPlot(app,inputsR,inputsL)
    startamp = 1;
    inputsR{6} = 1; %writeVTA
    inputsL{6} = 1; %writeVTA
    for i=1:9
        inputsR{2} = startamp;
        inputsL{2} = startamp;
        ea_generate_optim_vat(inputsR{:});
        ea_generate_optim_vat(inputsL{:});
        stimfolder = 'mA_33_14';
        Ihatvector(i,:) = predictImprovement(app,patient,stimfolder,1);
        amplitudevector(i) = startamp;
        startamp = startamp + 0.5;
    end
   figure;
   plot(amplitudevector,Ihatvector,'o-','linewidth',2,'markersize',5,'Color',[120/255,0,128/255]);

end

function options = setOPTS(patselect)
    options = ea_setopts_local;
    options.native = 0;
    options.groupmode = 1;
    options.groupid = 'cleartune';
    va = 0; % 0 for constant curr
    options = ea_getptopts(patselect, options);
    [coords_mm,trajectory,markers,elmodel,manually_corrected,coords_acpc]=ea_load_reconstruction(options);
    elstruct(1).coords_mm=coords_mm;
    elstruct(1).coords_acpc=coords_acpc;
    elstruct(1).trajectory=trajectory;
    elstruct(1).name = ['sub-', options.subj.subjId];
    elstruct(1).markers=markers;

    options.numcontacts=size(coords_mm{1},1);
    options.d3.verbose='off';
    options.d3.elrendering=1;	% hard code to viz electrodes in this setting.
    options.d3.exportBB=0;	% don't export brainbrowser struct by default
    options.d3.colorpointcloud=0;
    options.d3.hlactivecontacts=1;
    options.d3.showactivecontacts =1;
    options.d3.showpassivecontacts=1;
    options.d3.exportBB=0;
    options.expstatvat.do=0;
    options.leadprod = 'group';
    options.patient_list=patselect;
    options.d3.mirrorsides=0;
    options.atlasset = options.prefs.machine.vatsettings.horn_atlasset;
    options.patientname = ['sub-', options.subj.subjId];
return
end
function options=ea_setopts_local

    options.earoot=ea_getearoot;
    options.verbose=3;
    options.sides=1:2; % re-check this later..
    options.fiberthresh=1;
    options.writeoutstats=1;
    options.writeoutpm = 0;
    options.colormap='parula(64)';
    options.d3.write=1;
    options.d3.prolong_electrode=2;
    options.d3.writeatlases=1;
    options.macaquemodus=0;
return
end
