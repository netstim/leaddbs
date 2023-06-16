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
    reconstruction_file = dir([patlist{pt},filesep,'reconstruction',filesep,'*desc-reconstruction.mat']);
    reconstruction_file_path = [reconstruction_file.folder,filesep,reconstruction_file.name];
    [reconst, ~, ~, ~] = ea_get_reconstruction(reconstruction_file_path);
    % do not update S here, just get the bounds in mA!
    [min_bound_per_contact, max_bound_per_contact, ~] = ea_get_currents_per_contact(app.MinCylindricEditField_2.Value,app.MaxCylindricEditField.Value, app.MinSegmentedEditField.Value, app.MaxSegmentedEditField.Value, reconst, 0, 0);
    startptsR = zeros(1,app.inputVars.numContacts);
    startptsL = zeros(1,app.inputVars.numContacts);
    % set third contact (k2) to the middle of the higher current bound
    if abs(max_bound_per_contact(3)) > abs(min_bound_per_contact(3))
        startptsR(3) = max_bound_per_contact(3) / 2.0;
    else
        startptsR(3) = min_bound_per_contact(3) / 2.0;
    end
    startptsL(3) = startptsR(3);
    lbR = min_bound_per_contact;
    lbL = min_bound_per_contact;
    ubR = max_bound_per_contact;
    ubL = max_bound_per_contact;
    intergCond = zeros(1,app.inputVars.numContacts);
    A = repmat(-1,1,app.inputVars.numContacts);
    b = 5;
    Aeq = [];
    beq = [];

    optionsR=optimoptions('surrogateopt',...
        'ObjectiveLimit',-0.9,... % optimal solution with average Ihat ~0.9, lowest theoretical point is zero with an R of 1
        'MinSurrogatePoints',200,...
        'PlotFcn','surrogateoptplot',...
        'InitialPoints',startptsR,...
        'MaxFunctionEvaluations',1,...
        'Display','iter');
    optionsL=optimoptions('surrogateopt',...
        'ObjectiveLimit',-0.9,... % optimal solution with average Ihat ~0.9, lowest theoretical point is zero with an R of 1
        'MinSurrogatePoints',200,...
        'PlotFcn','surrogateoptplot',...
        'InitialPoints',startptsL,...
        'MaxFunctionEvaluations',1,...
        'Display','iter');

    %    'CheckpointFile',fullfile(fileparts(tractset.leadgroup),'optimize_status.mat'),...

    % check for parallel processing toolbox
    if ismember('Parallel Computing Toolbox',{toolboxes_installed.Name}) && useparallel
        optionsR.UseParallel=true;
        optionsL.UseParallel=true;
        parpool('Processes',2);
    end

    if all(ubR(:) <= 0.0)
        objconstrR=@(x)struct('Fval',nestedfunR_Monopolar(x,patlist{pt},ptindx));
        [XOptimR,fvalR,~,~,ipR]=surrogateopt(objconstrR,lbR,ubR,find(intergCond),A,b,Aeq,beq,optionsR);
    else
        objconstrR=@(x)struct('Fval',nestedfunR(x,patlist{pt},ptindx));
        [XOptimR,fvalR,~,~,ipR]=surrogateopt(objconstrR,lbR,ubR,optionsR);
    end

    writeVTA = 1;
    modelVTA = app.inputVars.modelVTA;
    %remove very small currents
    XOptimR(abs(XOptimR) < 0.000001) = 0.0;
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
    % new way to define inputs
    [inputsR,ampl_R,perc_val_R] = ea_get_inputs_for_optimizer(patlist{pt},XOptimR, modelVTA,writeVTA,1);
    ea_generate_optim_vat(inputsR{:});

    if all(ubL(:) <= 0.0)
        objconstrL=@(x)struct('Fval',nestedfunL_Monopolar(x,patlist{pt},ptindx));
        [XOptimL,fvalL,~,~,ipL]=surrogateopt(objconstrL,lbL,ubL,find(intergCond),A,b,Aeq,beq,optionsL);
    else
        objconstrL=@(x)struct('Fval',nestedfunL(x,patlist{pt},ptindx));
        [XOptimL,fvalL,~,~,ipL]=surrogateopt(objconstrL,lbL,ubL,optionsL);
    end

    % remove very small currents
    XOptimL(abs(XOptimL) < 0.000001) = 0.0;

    newoptimL = reformatX(XOptimL);
    if ~isempty(newoptimL)
        XOptimL = newoptimL;
    end
    writeVTA = 1;
    for i=1:size(ipL.X,1)
        newL = reformatX(ipL.X(i,:));
        if ~isempty(newL)
            ipL.X(i,:) = newL;
        end
    end
    save(fullfile(patlist{pt},'optimize_status_surrogate_L.mat'),'ipL','XOptimL');
    [inputsL,ampl_L,perc_val_L] = ea_get_inputs_for_optimizer(patlist{pt},XOptimL, app.inputVars.modelVTA,writeVTA,2);
    ea_generate_optim_vat(inputsL{:});
    options = setOPTS(patlist{pt});
    S = ea_initializeS(options);
    S = ea_cleartune_generateMfile([ampl_R,perc_val_R],[ampl_L,perc_val_L],S,0);
    save(fullfile(patlist{pt},'desc-stimparameters.mat'),'S');
    createIhatAmpPlot(app,inputsR,inputsL)
    avgIhat = ((-1*fvalR)+(-1*fvalL))/2;
    disp(['Optimal solution: Average Ihat(R,L) = ',num2str(avgIhat),'.']);
    warning on
    if ismember('Parallel Computing Toolbox',{toolboxes_installed.Name}) && useparallel
        poolobj = gcp('nocreate'); delete(poolobj);
    end

end
function [paramsR,paramsL] = genParams(app)

%for bipolar also add this
%bool_contact = repmat([0,1,1,0],app.inputVars.numContacts,1);
%bool_contact(app.inputVars.startContact,end) = 1;
%paramsR = [amplitude;perc_contactdist;bool_contact;case_contact];

perc_contactdist = repmat([-100, 0, 0, 0],app.inputVars.numContacts,1);
%case_contact = [100,100,1,100];
perc_contactdist(app.inputVars.startContact,end) = -100;

amplitude = [app.inputVars.minCurr, app.inputVars.maxCurr, 0, app.inputVars.minCurr];
paramsR = [amplitude;perc_contactdist];
paramsL = paramsR;

return
end

% Nested function that computes the objective function
function f = nestedfunR(X,pt,ptindx)

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
f.Fval=getFval(app,X,pt,1,ptindx);

return
end

function f = nestedfunL(X,pt,ptindx)
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
f.Fval=getFval(app,X,pt,2,ptindx);
% limit to 8mA

return
end


% Nested function that computes the objective function
function Fval = nestedfunR_Monopolar(X,pt,ptindx)

X = reformatX(X);

if isempty(X)
    Fval = Inf;
    return
end
tractsetclone=updateStim(tractsetclone,X);
Fval=getFval(app,X,pt,1,ptindx);

return
end

function Fval = nestedfunL_Monopolar(X,pt,ptindx)
X = reformatX(X);

if isempty(X)
    Fval = Inf;
    return
end

tractsetclone=updateStim(tractsetclone,X);
Fval=getFval(app,X,pt,2,ptindx);
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
% best fix ever
% for bipolar
%if the indicator is off, the contact should also be off
% for xx=10:17
%     if X(xx) == 0
%         X(xx-8) = 0;
%     end
% end
%doing it this way because harcoding indices is wrong and won't give us
%flexibility for future
%switch this on for bipolar setting
% whichpos = X > 0;
% whichpos(1) = 0;
% %whichpos(10:end) = 0;
% if any(whichpos)
%    sum_whichpos_r = sum(X(whichpos));
%    X(whichpos) = (X(whichpos)*100)/(sum_whichpos_r);
%    X(end) = 0;
% end
% else
%    X(end) = 100;
% end
% whichneg = X < 0;
% whichneg(1) = 0;
% % whichneg(10:end) = 0;
% if any(whichneg)
%     sum_whichneg_r = sum(X(whichneg));
%     X(whichneg) = (X(whichneg)*-100)/(sum_whichneg_r);
% end
% X(end) = X(end)*100;
% return
% end

function tractsetclone=updateStim(tractsetclone,X)
    disp('Parameters applied: ');
    %fprintf('Amplitude:%d\n',X(1));
    whichContact = find(X(1:end));
    fprintf('Active contact: k0%d\n',whichContact-1)
    fprintf('Current, mA: %d\n',X(1:end));
    %disp('Case activation:100%');
end

function Fval=getFval(app,X,pt,side,ptindx)
      
    %create a vta inside this function, send it to cleartune, get Ihat out and return it as
    %Fval
    
    % remove very small currents
    X(abs(X) < 0.000001) = 0.0;
    writeVTA = 0;

    % new way to define inputs
    inputs = ea_get_inputs_for_optimizer(pt,X, app.inputVars.modelVTA,writeVTA,side);
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
    [~,Ihat,actualimprovs] = runcrossval(app,'suggest',tractsetclone,{pt},side);
    %Ihat=[Ihat(1:length(Ihat)/2),Ihat(length(Ihat)/2+1:end)]; % Ihat is exported as a column vector for both sides. Reformat to Nx2.
    if sum(isnan(Ihat(:)))/length(Ihat(:))>0.3
        ea_warning(['Many (',num2str(sum(isnan(Ihat(:)))*100/length(Ihat(:))),' percent) stimulation settings were not covered well by the model. Stimulation suggestions may not be meaningful. Please adjust model parameters in ',...
            app.fibfiltmodelpath,'.']);
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
        Ihatvector(i) = predictImprovement(app,inputsR{2},stimfolder,1);
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
    options.patientname = options.subj.subjId;
return
end
function options=ea_setopts_local
    
    options.earoot=ea_getearoot;
    options.verbose=3;
    options.sides=1:2; % re-check this later..
    options.fiberthresh=1;
    options.writeoutstats=1;
    options.writeoutpm = 0;
    options.colormap=jet;
    options.d3.write=1;
    options.d3.prolong_electrode=2;
    options.d3.writeatlases=1;
    options.macaquemodus=0;
return
end
