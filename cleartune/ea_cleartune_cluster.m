function [tractsetclone]=ea_cleartune_cluster(tractsetclone,patlist,app,~,~)
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

for pt = 1:length(patlist)
    ptindx = pt;
    startptsR = app.startptsR;
    startptsL = app.startptsL;
    lbR = app.lbR; %min_bound_per_contactR;
    lbL = app.lbL; %min_bound_per_contactL;
    ubR = app.ubR; %max_bound_per_contactR;
    ubL = app.ubL; %max_bound_per_contactL;
    intergCond = zeros(1,app.inputVars.numContacts);
    A = repmat(-1,1,app.inputVars.numContacts);
    b = 5;
    Aeq = [];
    beq = [];

    reportpath = fullfile(patlist{pt},'log',filesep,'report.txt');
    fileID = fopen(reportpath,'a+');
    fprintf(fileID,['processing patient: ',patlist{pt}]);
    optionsR=optimoptions('surrogateopt',...
        'ObjectiveLimit',-0.9,... % optimal solution with average Ihat ~0.9, lowest theoretical point is zero with an R of 1
        'MinSurrogatePoints',200,...
        'PlotFcn','surrogateoptplot',...
        'InitialPoints',startptsR,...
        'MaxFunctionEvaluations',500,...
        'Display','iter');
    optionsL=optimoptions('surrogateopt',...
        'ObjectiveLimit',-0.9,... % optimal solution with average Ihat ~0.9, lowest theoretical point is zero with an R of 1
        'MinSurrogatePoints',200,...
        'PlotFcn','surrogateoptplot',...
        'InitialPoints',startptsL,...
        'MaxFunctionEvaluations',500,...
        'Display','iter');

    %    'CheckpointFile',fullfile(fileparts(tractset.leadgroup),'optimize_status.mat'),...

    % check for parallel processing toolbox
    if ismember('Parallel Computing Toolbox',{toolboxes_installed.Name}) && useparallel
        optionsR.UseParallel=true;
        optionsL.UseParallel=true;
        parpool('Processes',2);
    end

    % right side processing%
    disp('>>>>> Processing Right hemisphere..25% done >>>>>');
    fclose(fileID);
    if all(ubR(:) <= 0.0)
        objconstrR=@(x)struct('Fval',nestedfunR_Monopolar(x,patlist,ptindx,reportpath));
        [XOptimR,~,~,~,ipR]=surrogateopt(objconstrR,lbR,ubR,find(intergCond),A,b,Aeq,beq,optionsR);
    else
        objconstrR=@(x)struct('Fval',nestedfunR(x,patlist,ptindx,reportpath));
        [XOptimR,~,~,~,ipR]=surrogateopt(objconstrR,lbR,ubR,optionsR);
    end
    Flist = fopen('all');
    if ismember(fileID,Flist)
        fprintf(fileID,'\nRight side successfully completed!');
    else
        fileID = fopen(reportpath,'a+');
        fprintf(fileID,'\nRight side successfully completed!');
    end
    disp('>>>>> Right hemisphere completed..50% done >>>>>');

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
    fileID = fopen(reportpath,'a+');
    fprintf(fileID,'\nRight side optimizer paramters saved to patient directory!');
    fclose(fileID);
    [inputsR,ampl_R,perc_val_R] = ea_get_inputs_for_optimizer(patlist{pt},XOptimR, modelVTA,writeVTA,1);
    try
        ea_generate_optim_vat(inputsR{:});
        fprintf(fileID,'\nWriting out the winning VTA on the right side was successfull!');
        fprintf(fileID,'%d',inputsR{2});
        fclose(fileID);
    catch ME
        fileID = fopen(reportpath,'a+');
        fprintf(fileID,'\nWriting out the winning VTA on the right side did not succeed.\nHowever, your optimize status file is successfully saved.');
        fprintf(fileID,'%d',inputsR{2});
        fclose(fileID);
        disp(ME.message);
    end

   % left side processing%
    disp('>>>>> Processing Left hemisphere..75% done >>>>>');
    if all(ubL(:) <= 0.0)
        objconstrL=@(x)struct('Fval',nestedfunL_Monopolar(x,patlist,ptindx,reportpath));
        [XOptimL,~,~,~,ipL]=surrogateopt(objconstrL,lbL,ubL,find(intergCond),A,b,Aeq,beq,optionsL);
    else
        objconstrL=@(x)struct('Fval',nestedfunL(x,patlist,ptindx,reportpath));
        [XOptimL,~,~,~,ipL]=surrogateopt(objconstrL,lbL,ubL,optionsL);
    end
    Flist = fopen('all');
    if ismember(fileID,Flist)
        fprintf(fileID,'\nLeft side successfully completed!');
    else
        fileID = fopen(reportpath,'a+');
        fprintf(fileID,'\nLeft side successfully completed!');
    end
    % prune low currrent values
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
    fprintf(fileID,'\nLeft side optimizer paramters saved to patient directory!');
    [inputsL,ampl_L,perc_val_L] = ea_get_inputs_for_optimizer(patlist{pt},XOptimL, app.inputVars.modelVTA,writeVTA,2);
   try
       ea_generate_optim_vat(inputsL{:});
       fprintf(fileID,'\nWriting out the winning VTA on the left side was successfull!');
       fprintf(fileID,'%d',inputsL{2});
   catch ME
        disp(ME.message);
        fprintf(fileID,'\nWriting out the winning VTA on the left side did not succeed.\nHowever, your optimize status file is successfully saved.');
        fprintf(fileID,'%d',inputsL{2});
   end
   disp('>>>>> Left hemisphere completed..100% done >>>>>');
   
   %finally save the stimulation parameters in the patient directory
   options = setOPTS(patlist{pt});
   S = ea_initializeS(options);
   S = ea_cleartune_generateMfile([ampl_R,perc_val_R],[ampl_L,perc_val_L],S,0);
   save(fullfile(patlist{pt},'desc-stimparameters.mat'),'S');
   disp('>>>>> Stimulation parameters saved successfully! >>>>>');
   fprintf(fileID,'\nStimulation parameters saved!');
   fclose(fileID);
    %finally, switch off parallel processing
    warning on
    if ismember('Parallel Computing Toolbox',{toolboxes_installed.Name}) && useparallel
        poolobj = gcp('nocreate'); delete(poolobj);
    end
    disp("Process done ***")

end


% Nested function that computes the objective function
function Fval = nestedfunR(X,patlist,ptindx,reportpath)
side = 1; %right side
% limit to 8mA
if sum(X(X>0)) > sum(abs(X(X<0)))
    Ineq = sum(X(X>0)) - 5.0;
else
    Ineq = sum(abs(X(X<0))) - 5.0;
end

X = reformatX(X);

if isempty(X)
    Fval = Inf;
    return
end
tractsetclone=updateStim(tractsetclone,X,side);
[Fval]=getFval(app,X,patlist,side,ptindx,reportpath);
return
end

function Fval = nestedfunL(X,patlist,ptindx,fileID)
side = 2; %left side
X = reformatX(X);

if sum(X(X>0)) > sum(abs(X(X<0)))
    Ineq = sum(X(X>0)) - 8.0;
else
    Ineq = sum(abs(X(X<0))) - 8.0;
end

if isempty(X) 
    Fval = Inf;
    return
end
tractsetclone=updateStim(tractsetclone,X,side);
Fval=getFval(app,X,patlist,side,ptindx,fileID);
% limit to 8mA
return
end

         
% Nested function that computes the objective function
function Fval = nestedfunR_Monopolar(X,patlist,ptindx,fileID)
side = 1; %side = 1 since it is right side
X = reformatX(X);

if isempty(X)
    Fval = Inf;
    fprintf(fileID,'\nEmpty input variables. Returning Fval of Inf.');
    return
end
tractsetclone=updateStim(tractsetclone,X,side);
Fval=getFval(app,X,patlist,side,ptindx,reportpath);
return
end

function Fval = nestedfunL_Monopolar(X,patlist,ptindx,reportpath)
side = 2; %since its left side
X = reformatX(X);

if isempty(X)
    Fval = Inf;
    return
end

tractsetclone=updateStim(tractsetclone,X,side);
Fval=getFval(app,X,patlist,side,ptindx,reportpath);

% limit to 8mA
return
end


function X = reformatX(ipX)
if all(~ipX) %all contacts have zero % activation its not allowed
    X = [];
    return
end


thresh = abs(app.inputVars.MinCylindricCurr/10);
%normalization, which will ruin the function progression
%whichneg = abs(X) < (abs(app.MinCylindricEditField_2.Value/100) &&;
for yy = 1:length(ipX)
    if ipX(yy) ~= 0 && abs(ipX(yy)) < thresh
        whichneg(yy) = 1;
    else
        whichneg(yy) = 0;
    end
end
sum_X = sum(ipX);
if any(whichneg)
    ipX(find(whichneg)) = 0;
    sum_newX = sum(ipX);
    normalizer = (sum_X/sum_newX);
    X = (ipX.*normalizer);
else
    X = ipX;
end
return
end

function tractsetclone=updateStim(tractsetclone,X,side)
        disp('Parameters applied: ');
        whichContact = find(X(1:end));
        if side == 1
            fprintf('Active contact: k0%d\n',whichContact-1);
        else
            if whichContact < 2
                fprintf('Active contact: k0%d\n',whichContact+8);
            else
                fprintf('Active contact: k%d\n',whichContact+8);
            end
        end
        fprintf('Current, mA: %d\n',X(1:end));
end

function Fval=getFval(app,X,patlist,side,ptindx,reportpath)
        %create a vta inside this function, send it to cleartune, get Ihat out and return it as
        %Fval
        fileID = fopen(reportpath,'a+');
        selpat = patlist{ptindx};
        % remove very small currents
        X(abs(X) < 0.000001) = 0.0;
        writeVTA = 0;

        % new way to define inputs
        inputs = ea_get_inputs_for_optimizer(selpat,X, app.inputVars.modelVTA,writeVTA,side);
        try
            [Efields,allS]=ea_generate_optim_vat(inputs{:});
        catch ME
            disp(ME.message);
            Fval = NaN;
            fprintf(fileID,'\nComputation of VTA failed');
            fprintf(fileID,'%d',[inputs{2} ', side:' side]);
            return
        end
        app.protocol{ptindx}.inputs=inputs;
        app.protocol{ptindx}.Efields=Efields;
        app.protocol{ptindx}.allS=allS;

        try
            [~,Ihat,actualimprovs] = runcrossval(app,'suggest',tractsetclone,patlist,ptindx,side);
            if any(isnan(Ihat))
                nanidx = find(isnan(Ihat));
                fprintf(fileID,'\n%s%d',['Calculation of Ihat returned NaN for the following idx: ',nanidx]);
            end
            preFval = calculateFval(app,Ihat,actualimprovs,side,ptindx);
            Fval = -1*preFval;
            %add penalty function if user chooses
            Fval = penaltyFunc(app,X,Fval);
        catch
            fprintf(fileID,'\n%s','Calculation of Ihat failed for patient:');
            fprintf(fileID,'\n%s',patlist{ptindx});
            fprintf(fileID,'\n%s%d',['side: ',side]);
            disp("Program exited with error: Ihat not calculated");
        end
        try
            fclose(fileID);
        end
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
    if strcmp(app.Applypenaltyusing,'Quadratic curve')
        xx = app.PenaltyVal;
        yy = abs(sum(X));
        Fval = Fval + ((yy-xx)^2).*0.02;
    elseif strcmp(app.Applypenaltyusing,'Range of amplitudes')
        if sum(abs(X)) > 3 || sum(abs(X)) < 2
            Fval = Fval + app.ApplyapenaltyvalueofEditField.Value;
        end
    else
        return 
    end
    return
end
function createIhatAmpPlot(app,inputsR,inputsL)
    startamp = 1;
    inputsR{6} = 1; %writeVTA
    inputsL{6} = 1; %writeVTA
    patient = {inputsR{1}};

    for i=1:9
        inputsR{2} = startamp;
        inputsL{2} = startamp;
        ea_generate_optim_vat(inputsR{:});
        ea_generate_optim_vat(inputsL{:});
        stimfolderR = 'mA_33_14_side_1';
        Ihatvector(i,1) = predictImprovement(app,patient,stimfolderR,1);
        stimfolderL = 'mA_33_14_side_2';
        Ihatvector(i,2) = predictImprovement(app,patient,stimfolderL,1);
        amplitudevector(i) = startamp;
        startamp = startamp + 0.5;
    end
   figure;
   plot(amplitudevector,Ihatvector(:,1),'o-','linewidth',2,'markersize',5,'Color',[120/255,0,128/255]);
   hold on
   plot(amplitudevector,Ihatvector(:,2),'o-','linewidth',2,'markersize',5,'Color',[135/255,206/255,200/255]);
   filename = fullfile(inputsR{1},'Ihat_vs_amplitude.png');
   saveas(figure,filename);

   
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

function [I,Ihat,actualimprovs] = runcrossval(app,event,tractsetclone,patlist,ptidx,side)
tractsetclone.calculate_cleartune(app.protocol{ptidx}.Efields);
app.tractset.cleartuneresults=tractsetclone.cleartuneresults; % this is fine to copy to the main object for visualization later.

%% >>>>>>>>>> BEGIN INJECTION:

%% temporarily append cleartune results to results to feed into
% crossvalidation:

% this is crucial to be cleaned up again and cannot go
% wrong, i.e. anything in between injection and cleanup
% needs to be in a try/catch statement.

if ~isempty(tractsetclone.cleartuneinjected)
    if ~strcmp(tractsetclone.cleartuneinjected.status,'clean')
        ea_error('Fiber Filtering File is corrupted. Aborting.');
    end
end

tractsetclone.cleartuneinjected.status='injected';
gpatsel=tractsetclone.patientselection;
if tractsetclone.mirrorsides
    gpatsel=[gpatsel,gpatsel + length(tractsetclone.M.patient.list)]; %this should not be changed.
end
%find length of selected patients to add the new guy in
origlen=size(tractsetclone.results.(ea_conn2connid(tractsetclone.connectome)).(ea_method2methodid(tractsetclone)).fibsval{1}(:,gpatsel),2)/2;
%to find index of the new guy from cleartune.results
numEfields = size(app.protocol{ptidx}.Efields(side),1);
origix=1:origlen; %the last value of this array will indicate where the new efield will be added.
efieldix=1:numEfields;
%first make a copy of the results
tmpfibvals=tractsetclone.results.(ea_conn2connid(tractsetclone.connectome)).(ea_method2methodid(tractsetclone)).fibsval;
%find the last patient to insert new guy into
lastpt = origix(end);
%this will also be needed for array indexing
allpatslen = length(tractsetclone.M.patient.list);
for hem=1:2
    tractsetclone.results.(ea_conn2connid(tractsetclone.connectome)).(ea_method2methodid(tractsetclone)).fibsval{hem}=...
        [tractsetclone.results.(ea_conn2connid(tractsetclone.connectome)).(ea_method2methodid(tractsetclone)).fibsval{hem}(:,1:lastpt),...
        tractsetclone.cleartuneresults.(ea_conn2connid(tractsetclone.connectome)).(ea_method2methodid(tractsetclone)).fibsval{hem}(:,efieldix),...
        tractsetclone.results.(ea_conn2connid(tractsetclone.connectome)).(ea_method2methodid(tractsetclone)).fibsval{hem}(:,lastpt+1:lastpt+allpatslen),...
        tractsetclone.cleartuneresults.(ea_conn2connid(tractsetclone.connectome)).(ea_method2methodid(tractsetclone)).fibsval{hem}(:,efieldix),...
        tractsetclone.results.(ea_conn2connid(tractsetclone.connectome)).(ea_method2methodid(tractsetclone)).fibsval{hem}(:,lastpt+allpatslen+1:end)];
end
tractsetclone.allpatients=[tractsetclone.M.patient.list;repmat({''},size(app.protocol{ptidx}.Efields,1),1)];

tractsetclone.responsevar=[tractsetclone.responsevar;nan(size(app.protocol{ptidx}.Efields,1),size(tractsetclone.responsevar,2))];
reduceresponsevar=0;
switch event
    case 'suggest'
        if size(tractsetclone.responsevar,2) == 1 % global scores
            tractsetclone.responsevar=[tractsetclone.responsevar,tractsetclone.responsevar];    % we will need maximal stim settings for each site separately, so we will pretend responsevar is bihemispheric
            reduceresponsevar=1;
        end
end
switch tractsetclone.multitractmode
    case 'Single Tract Analysis'
    otherwise
        for wv=1:length(tractsetclone.subscore.vars)
            tractsetclone.subscore.weightvars{wv} = [tractsetclone.subscore.weightvars{wv}(:);repmat(app.symptomWeightVar{ptidx,side}(wv),size(app.protocol{ptidx}.Efields,1),1)];
            tractsetclone.subscore.vars{wv}=[tractsetclone.subscore.vars{wv}(:);ones(size(app.protocol{ptidx}.Efields,1),1)]; % need to add in weights here.
        end

        app.trainings=repmat({app.tractset.patientselection},length(patlist),1);

        try % this has to be in a try/catch because we still have to cleanup the fibfilt again if something would go wrong (next section).
            cvp.NumTestSets = 1;
            % training and test indices from the items list
            test = origlen+1:origlen+size(app.protocol{ptidx}.Efields,1);
            % Patient selected based on the training and test indices
            tractsetclone.customselection = unique([app.trainings{ptidx}, test]);
            cvp.training{1} = ismember(tractsetclone.customselection, app.trainings{ptidx});
            cvp.test{1} = ismember(tractsetclone.customselection, test);
            % Construct cvp struct
            [I, Ihat, actualimprovs]=tractsetclone.crossval(cvp); % actualimprovs is for each side separately.

            %remove results again
            for hem=1:2
                tractsetclone.results.(ea_conn2connid(tractsetclone.connectome)).(ea_method2methodid(tractsetclone)).fibsval{hem}=...
                    [tractsetclone.results.(ea_conn2connid(tractsetclone.connectome)).(ea_method2methodid(tractsetclone)).fibsval{hem}(:,1:lastpt),...
                    tractsetclone.results.(ea_conn2connid(tractsetclone.connectome)).(ea_method2methodid(tractsetclone)).fibsval{hem}(:,lastpt+2:lastpt+allpatslen+1),... %remove efield of new guy
                    tractsetclone.results.(ea_conn2connid(tractsetclone.connectome)).(ea_method2methodid(tractsetclone)).fibsval{hem}(:,lastpt+allpatslen+3:end)]; %remove mirrorred of new guy
            end
            if isequaln(tractsetclone.results.(ea_conn2connid(tractsetclone.connectome)).(ea_method2methodid(tractsetclone)).fibsval,tmpfibvals)
                tractsetclone.cleartuneinjected.status='clean';
            else
                ea_error('Something went wrong.');
            end
            return
        end
end
end

