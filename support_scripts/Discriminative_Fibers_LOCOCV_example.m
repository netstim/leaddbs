load('/path/to/your/LEAD_groupanalysis.mat'); % this creates variable M
M.ui.connectomename='HCP_MGH_30fold_groupconnectome (Horn 2017)'; % ... which can be modified afterwards
M.ui.listselect=1:length(M.patient.list); % usually helpful to create the model for all patients but then use only the ones defined in the script to run the analyses.
prefs=ea_prefs;

%% static part - usually no need to edit - can edit it by configuring lead group file correctly and closing it.
discfiberssetting = prefs.machine.lg.discfibers;
[fibsweighted,fibsin,fibsval,iaix]=ea_discfibers_calcdiscfibers(M,discfiberssetting);


%% flexible part - how to set up prediction - example is leave one cohort out crossvalidation:
allpts=1:length(M.patient.list);
I=M.clinical.vars{M.ui.clinicallist};
ea_dispercent(0,'Predicting left-out patients');
nfibsval=fibsval; nfibsval(nfibsval==0)=nan; % only used in spearmans correlations

%nfibsval(nfibsval<50)=nan; % could be used to set a threshold on the E-Field method.
if discfiberssetting.statmetric==1
    % In t-test method, make sure fibers are at least connected to 0.2 percent of VTAs and not
    % more than to 0.8 percent.
    discthresh=0.2; % could be changed but should not exceed ~0.45
    fibsval(sum(fibsval,2)<round(discthresh*size(fibsval,2)),:)=0;
    fibsval(sum(fibsval,2)>round((1-discthresh)*size(fibsval,2)),:)=0;
end

for group=unique(M.patient.group)'
    opts=allpts; opts(M.patient.group(allpts)==group)=[]; % generate variable of patients on which model will be built.
    switch discfiberssetting.statmetric
        case 1 % ttests, vtas - see Baldermann et al. 2019 Biological Psychiatry
            optsval=fibsval(:,opts); % all other patients connections to each fibertract
            allvals=repmat(I(opts)',size(optsval,1),1); % improvement values (taken from Lead group file or specified in line 12).
            fibsimpval=allvals; % Make a copy to denote improvements of connected fibers
            fibsimpval(~logical(optsval))=nan; % Delete all unconnected values
            nfibsimpval=allvals; % Make a copy to denote improvements of unconnected fibers
            nfibsimpval(logical(optsval))=nan; % Delete all connected values
            [~,~,~,Model]=ttest2(fibsimpval',nfibsimpval'); % Run two-sample t-test across connected / unconnected values
            Model.tstat(p>0.5)=nan; % discard noisy fibers (optional or could be adapted)
            for pt=find(M.patient.group==group)
                thisptval=fibsval(:,pt); % this patients connections to each fibertract (1 = connected, 0 = unconnected) 
                Ihat(pt)=ea_nansum(Model.tstat'.*thisptval); % I hat is the estimate of improvements (not scaled to real improvements)
            end
        case 2 % spearmans correlations, efields - see Li et al. 2019 TBP
            Model=corr(nfibsval(:,opts)',I(opts),'rows','pairwise','type','Spearman'); % generate optimality values on all but left out patients
            for pt=find(M.patient.group(allpts)==group)
                Ihat(pt)=ea_nansum(Model.*nfibsval(:,pt)); % I hat is the estimate of improvements (not scaled to real improvements)
            end
    end
    ea_dispercent(pt/length(allpts));
end
ea_dispercent(1,'end');
loginx=zeros(size(Ihat)); loginx(allpts)=1;
Ihat(~loginx)=nan; % make sure info of not included patients are not used


h=ea_corrplot(I,Ihat',{'Disc. Fiber prediction LOOCV','Empirical','Predicted'},'permutation_spearman');
saveas(h,'my_result.png');






